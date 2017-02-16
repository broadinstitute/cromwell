package cromwell.core.callcaching.docker.registryv2

import akka.NotUsed
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.HttpHeader.ParsingResult
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge}
import cromwell.core.callcaching.docker.DockerHashActor._
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow._
import cromwell.core.callcaching.docker.registryv2.flows.{FlowUtils, HttpFlowWithRetry}
import cromwell.core.callcaching.docker.{DockerFlow, DockerHashResult, DockerImageIdentifierWithoutHash}
import spray.json._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object DockerRegistryV2AbstractFlow {
  type HttpDockerFlow = Flow[(HttpRequest, (DockerHashContext, HttpRequest)), (Try[HttpResponse], (DockerHashContext, HttpRequest)), NotUsed]
  
  val DigestHeaderName = "Docker-Content-Digest".toLowerCase
  val AcceptHeader = HttpHeader.parse("Accept", "application/vnd.docker.distribution.manifest.v2+json") match {
    case ParsingResult.Ok(header, errors) if errors.isEmpty => header
    case ParsingResult.Ok(header, errors) =>
      errors foreach { err => logger.warn(err.formatPretty) }
      header
    case ParsingResult.Error(error) => throw new RuntimeException(error.formatPretty)
  }
}

/**
  * Implements logic to build a flow that can retrieve a docker hash 
  * from a registry server conforming to the docker registry API v2 Specification:
  * https://docs.docker.com/registry/spec/api/
  */
abstract class DockerRegistryV2AbstractFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer) extends DockerFlow {
  // Wraps the Http flow in a retryable flow to enable auto retries
  final private val httpFlowWithRetry = new HttpFlowWithRetry[DockerHashContext](httpClientFlow, isRetryable = isRetryable).flow

  final private val tokenFlow = {
    val responseHandlerFlow = Flow[(HttpResponse, DockerHashContext)].mapAsync(1)(Function.tupled(tokenResponseHandler))
    requestTransformFlow(buildTokenRequest, responseHandlerFlow)
  }

  final private val manifestFlow = {
    val responseHandlerFlow = Flow[(HttpResponse, DockerHashContext)].map(Function.tupled(manifestResponseHandler))
    requestTransformFlow(Function.tupled(buildManifestRequest _), responseHandlerFlow)
  }
  
  /**
    * Builds a flow for this docker registry
    */
  def buildFlow() = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._
    
    // Decouple the token flow from the manifest flow with ".async" 
    // this way they can run in parallel from each other 
    // while still maintaining the final ordering
    val token = builder.add(tokenFlow.async)
    val manifest = builder.add(manifestFlow.async)
    val mergeResponses = builder.add(Merge[(DockerHashResponse, DockerHashContext)](3))
    
    token.out0 ~> manifest.in
                  manifest.out0 ~> mergeResponses.in(0)
                  manifest.out1 ~> mergeResponses.in(1)
    token.out1       ~>            mergeResponses.in(2)
    
    FlowShape(token.in, mergeResponses.out)
  }

  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = dockerImageIdentifierWithoutHash.host.contains(registryHostName)
  
  /* Methods that must to be implemented by a subclass */

  /**
    * (e.g registry-1.docker.io)
    */
  protected def registryHostName: String

  /**
    * (e.g auth.docker.io)
    */
  protected def authorizationServerHostName: String

  /**
    * Builds the list of headers for the token request
    */
  protected def buildTokenRequestHeaders(dockerHashContext: DockerHashContext): List[HttpHeader]
  
  /* Methods that may be overridden by a subclass */

  /**
    * Returns true if an http response failed in such a way that it can be retried,
    * false otherwise.
    */
  protected def isRetryable(httpResponse: HttpResponse) = {
    httpResponse.status match {
      case StatusCodes.TooManyRequests => true
      case StatusCodes.InternalServerError => true
      case _ => false
    }
  }

  /**
    * service parameter that can be used in the token request
    * (e.g https://auth.docker.io/token?service=registry.docker.io&scope=repository:library/ubuntu:pull)
    */
  protected def serviceName: Option[String] = None

  /**
    * Builds the token URI to be queried based on a DockerImageID
    */
  protected def buildTokenRequestUri(dockerImageID: DockerImageIdentifierWithoutHash): String = {
    val service = serviceName map { name => s"service=$name&" } getOrElse ""
    s"https://$authorizationServerHostName/token?${service}scope=repository:${dockerImageID.name}:pull"
  }

  /**
    * Http method used for the manifest request
    */
  protected def manifestRequestHttpMethod: HttpMethod = HttpMethods.GET

  /**
    * Maps a failed http response to a DockerHashFailureResponse
    */
  protected def httpFailureMapper(response: HttpResponse, dockerHashContext: DockerHashContext) = {
    response.status match {
      case StatusCodes.NotFound => (DockerHashNotFound(dockerHashContext.dockerImageID), dockerHashContext)
      case StatusCodes.Unauthorized => (DockerHashUnauthorized(dockerHashContext.dockerImageID), dockerHashContext)
      case other => (
        DockerHashFailedResponse(
          new Exception(s"Docker hash lookup failed with code ${other.intValue()}. ${other.defaultMessage()}"), dockerHashContext.dockerImageID
        ),
        dockerHashContext
        )
    }
  }

  /**
    * Generic method to build a flow that creates a request, sends it,
    * gets the result and performs some transformation on the result
    * @param requestBuilderFunction Function that builds a request from an input object of type A
    * @param responseHandlerFlow Function that transforms a response in an object of type Try[B]
    * @tparam A type of the input value of the flow
    * @tparam B type of the output value of the flow, if successful
    * @return flow
    */
  private def requestTransformFlow[A, B](
                                          requestBuilderFunction: A => (HttpRequest, DockerHashContext),
                                          responseHandlerFlow: Flow[(HttpResponse, DockerHashContext), (Try[B], DockerHashContext), NotUsed]
                                        ) = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // Build a request from an A
    val request = builder.add(Flow.fromFunction(requestBuilderFunction))

    // Perform http requests
    val http = builder.add(httpFlowWithRetry)
    val successfulResponse = http.out0.via(responseHandlerFlow)
    val failedResponse = http.out1.map(Function.tupled(httpFailureMapper))

    // Partition parsing results
    val parser = builder.add(FlowUtils.fanOutTry[B, DockerHashContext])
    val successfulParsing: Outlet[(B, DockerHashContext)] = parser.out0
    val failedParsing = parser.out1.map(Function.tupled(throwableToFailureMessage))

    // Merge failures from failed response and failed parsing
    val mergeFailures = builder.add(Merge[(DockerHashResponse, DockerHashContext)](2))

    request ~> http.in
    successfulResponse ~> parser.in

    failedResponse ~> mergeFailures.in(0)
    failedParsing ~> mergeFailures.in(1)

    new FanOutShape2(request.in, successfulParsing, mergeFailures.out)
  }

  /**
    * Builds the token request
    */
  private def buildTokenRequest(dockerHashContext: DockerHashContext) = {
    val authorizationHeaders = buildTokenRequestHeaders(dockerHashContext)

    val tokenHttpRequest = HttpRequest(
      method = HttpMethods.GET,
      uri = buildTokenRequestUri(dockerHashContext.dockerImageID),
      headers = authorizationHeaders
    )

    (tokenHttpRequest, dockerHashContext)
  }

  /**
    * Parse the http response coming back from the token request and returns the body as JsObject
    */
  private def tokenResponseHandler(response: HttpResponse, dockerHashContext: DockerHashContext) = response match {
    case httpResponse if httpResponse.status.isSuccess() =>
      Unmarshal(httpResponse.entity).to[String] map { body =>
        (
          Try(body.parseJson.asJsObject) flatMap extractToken,
          dockerHashContext
        )
      } recoverWith {
        case failure => Future.successful((Failure(failure), dockerHashContext))
      }
    case httpResponse =>
      Future.successful(
        (
        Failure(new Exception(s"Request to obtain docker registry access token returned with a failure code: ${httpResponse.status.value}")),
        dockerHashContext
        )
      )
  }

  /**
    * Extract the access token from the json body of the http response
    */
  private def extractToken(jsObject: JsObject) = {
    jsObject.fields.get("token") match {
      case Some(token: JsString) => Success(token.value)
      case Some(other) => Failure(new Exception("Token response contains a non-string token field"))
      case None => Failure(new Exception("Token response did not contain a token field"))
    }
  }

  /**
    * Builds the manifest URI to be queried based on a DockerImageID
    */
  private def buildManifestUri(dockerImageID: DockerImageIdentifierWithoutHash): String = {
    s"https://$registryHostName/v2/${dockerImageID.name}/manifests/${dockerImageID.reference}"
  }

  /**
    * Builds the manifest http request
    */
  private def buildManifestRequest(token: String, dockerHashContext: DockerHashContext) = {
    val authorizationHeader = Authorization(OAuth2BearerToken(token))

    val manifestRequest = HttpRequest(
      method = manifestRequestHttpMethod,
      uri = buildManifestUri(dockerHashContext.dockerImageID),
      headers = scala.collection.immutable.Seq(AcceptHeader, authorizationHeader)
    )

    (manifestRequest, dockerHashContext)
  }

  /**
    * Handles the http response from the manifest request
    */
  private def manifestResponseHandler(response: HttpResponse, dockerHashContext: DockerHashContext) = {
    (extractContentDigest(response.headers), dockerHashContext) 
  }
  
  /**
    * Extracts the digest from the response headers
    */
  private def extractContentDigest(headers: Seq[HttpHeader]) = {
    headers find { _.is(DigestHeaderName) } match {
      case Some(digestHeader) => Success(DockerHashResponseSuccess(DockerHashResult(digestHeader.value())))
      case None => Failure(new Exception("Cannot find digest header"))
    }
  }

  /**
    * Wraps a throwable into a DockerHashFailedResponse
    */
  private def throwableToFailureMessage(throwable: Throwable, dockerHashContext: DockerHashContext): (DockerHashResponse, DockerHashContext) = {
    (DockerHashFailedResponse(throwable, dockerHashContext.dockerImageID), dockerHashContext)
  }
}
