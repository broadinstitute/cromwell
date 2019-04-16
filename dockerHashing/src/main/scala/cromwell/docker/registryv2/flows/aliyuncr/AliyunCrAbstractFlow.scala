package cromwell.docker.registryv2.flows.aliyuncr

import akka.NotUsed
//import akka.actor.Scheduler
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
// import akka.http.scaladsl.model.HttpHeader.ParsingResult
import akka.http.scaladsl.model._
// import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge}
import akka.stream.FlowShape
import cromwell.docker.DockerHashActor._
import cromwell.docker.registryv2.flows.aliyuncr.AliyunCrAbstractFlow._
import cromwell.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest
import cromwell.docker.registryv2.flows.{FlowUtils}
import cromwell.docker.{DockerFlow, DockerHashResult, DockerImageIdentifierWithoutHash}
import spray.json._
import DefaultJsonProtocol._

//import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
//import scala.language.postfixOps
import scala.util.{Failure, Try}
import scala.util.matching.Regex

import com.aliyuncs.auth.{BasicCredentials}
import com.aliyuncs.DefaultAcsClient;
import com.aliyuncs.cr.model.v20160607.GetRepoTagsRequest;
// import com.aliyuncs.exceptions.ClientException;
import com.aliyuncs.http.{HttpResponse => JavaHttpResponse};
import com.aliyuncs.profile.DefaultProfile;
import com.aliyuncs.profile.IClientProfile;

object AliyunCrAbstractFlow {
  type HttpDockerFlow = Flow[(HttpRequest, ContextWithRequest[DockerHashContext]), (Try[HttpResponse], ContextWithRequest[DockerHashContext]), NotUsed]
  //val StrictTimeout = 30 seconds

  val ProductName = "cr"

  val HashAlg = "sha256"
  
  //val DigestHeaderName = "Docker-Content-Digest".toLowerCase

  //response-body-data-manifest-fsLayers
  // val AcceptHeader = HttpHeader.parse("Accept", "application/vnd.docker.distribution.manifest.v2+json") match {
  //   case ParsingResult.Ok(header, errors) if errors.isEmpty => header
  //   case ParsingResult.Ok(header, errors) =>
  //     errors foreach { err => logger.warn(err.formatPretty) }
  //     header
  //   case ParsingResult.Error(error) => throw new RuntimeException(error.formatPretty)
  // }

  case class UnsuccessfulHttpResponseException(httpResponse: HttpResponse) extends Exception
  case class UnsuccessfulAliyunHttpResponseException(httpResponse: JavaHttpResponse) extends Exception

}

/**
  * Implements logic to build a flow that can retrieve a docker hash 
  * from a Aliyun CR
  */
abstract class AliyunCrAbstractFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer/*, notUsed: Scheduler*/) extends DockerFlow {
  // Wraps the Http flow in a retryable flow to enable auto retries
  // final private val httpFlowWithRetry = GraphDSL.create() { implicit builder =>
  //   import GraphDSL.Implicits._
  //   val retryHttpFlow = builder.add(new HttpFlowWithRetry[DockerHashContext](httpClientFlow).flow)
    
  //   // Force the response entity to be strict. This makes sure that whatever happens to the response later we
  //   // won't leave it hanging and potentially lock the pool. 
  //   // See http://doc.akka.io/docs/akka-http/10.0.5/scala/http/client-side/request-level.html#using-the-future-based-api-in-actors
  //   // Note that in this particular case it's ok to force the loading of the entity in memory
  //   // because we use HEAD Http method when we only care about the headers. Therefore there's no unnecessary memory usage.
  //   /* Returns a (Try[HttpResponse], DockerHashContext) */
  //   val strictHttpResponse = retryHttpFlow.out0.mapAsync(1){
  //     case (response, context) => response.toStrict(StrictTimeout) map { Success(_) -> context } recoverWith {
  //       case failure => Future.successful(Failure(failure) -> context)
  //     }
  //   }
    
  //   // Splits successful `toStrict` responses from failures
  //   val partitionStrictResponse = builder.add(FlowUtils.fanOutTry[HttpResponse, DockerHashContext])

  //   // Merge failures from retryHttpFlow.out1 (failed http responses)
  //   // and partitionStrictResponse.out1 (failed to `toStrict` the response)
  //   val mergeFailures = builder.add(Merge[(Throwable, DockerHashContext)](2))
    
  //   strictHttpResponse.outlet ~> partitionStrictResponse.in

  //   retryHttpFlow.out1 ~> mergeFailures
  //   partitionStrictResponse.out1 ~> mergeFailures
    
  //   new FanOutShape2(retryHttpFlow.in, partitionStrictResponse.out0, mergeFailures.out)
  // }

  // private [registryv2] val tokenFlow = {
  //   val responseHandlerFlow = Flow[(HttpResponse, DockerHashContext)]
  //     .mapAsync(1)(Function.tupled(tokenResponseHandler))
  //     // Map the Try[String] token to a Try[Option[String]].
  //     // This allows potential users of this class to override this flow and return an empty token
  //     .map { case (tryToken, context) => (tryToken map Option.apply, context) }
    
  //   requestTransformFlow(buildTokenRequest, responseHandlerFlow)
  // }

  // private [registryv2] val manifestFlow = GraphDSL.create() { implicit builder =>
  //   val maniFlow = builder.add(Flow.fromFunction(Function.tupled(getManifest _)))
  //   val emptyFailurePort = builder.add(Source.empty[(DockerHashResponse, DockerHashContext)])
  //   new FanOutShape2(maniFlow.in, maniFlow.out, emptyFailurePort.out)
  // }
  
  /**
    * Builds a flow for this docker registry
    */
  def buildFlow() = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // val token = builder.add(tokenFlow)
    // // val manifest = builder.add(Flow.fromFunction(Function.tupled(getManifest _)))
    // // val mergeResponses = builder.add(Merge[(DockerHashResponse, DockerHashContext)](2))
    
    // // token.out0 ~> manifest ~> mergeResponses.in(0)
    // // token.out1       ~>       mergeResponses.in(1)
    // val manifest = builder.add(manifestFlow)
    // val mergeResponses = builder.add(Merge[(DockerHashResponse, DockerHashContext)](3))
    
    // token.out0 ~> manifest.in
    //               manifest.out0 ~> mergeResponses.in(0)
    //               manifest.out1 ~> mergeResponses.in(1)
    // token.out1       ~>            mergeResponses.in(2)
    // FlowShape(token.in, mergeResponses.out)

    //  val responseHandlerFlow = Flow[(HttpResponse, DockerHashContext)]
    //   .mapAsync(1)(Function.tupled(tokenResponseHandler))
    //   .map { case (tryToken, context) => (tryToken map Option.apply, context) }

    //  val manifestFlow = Flow[(DockerHashSuccessResponse, DockerHashContext)]
    //   .mapAsync(1)((getManifest _))
    //   .map { case (hashResponse, context) => (hashResponse map Option.apply, context) }

    // val manifest = builder.add(manifestFlow)

    val manifest = builder.add(
      Flow
        .fromFunction(getManifest _)
        .mapAsync(1)(identity)
    )

    //     // Check if the hash was not found
    // val partition = builder.add(Partition[(Try[DockerHashResponse], DockerHashContext)](2, {
    //   case (Success(_: DockerHashSuccessResponse), _) => 0 // If we didn't find the docker hash, go to out(0)
    //   case _ => 1 // Otherwise, go to out(1)
    // }))

    // Partition processed responses
    val processResponse = builder.add(FlowUtils.fanOutTry[DockerHashResponse, DockerHashContext])
    
    // val successfulProcessing = processResponse.out0
    // val failedProcessing = processResponse.out1.map(Function.tupled(throwableToFailureMessage))

    val merge = builder.add(Merge[(DockerHashResponse, DockerHashContext)](2))



    manifest.out ~> processResponse.in
                    processResponse.out0                           ~>                        merge.in(0)
                    processResponse.out1.map(Function.tupled(throwableToFailureMessage)) ~>  merge.in(1)
    
    FlowShape(manifest.in, merge.out)
  }

  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash) = false
  
  // /* Methods that must to be implemented by a subclass */

  // /**
  //   * (e.g registry-1.docker.io)
  //   */
  // protected def registryHostName: String

  // /**
  //   * (e.g auth.docker.io)
  //   */
  // protected def authorizationServerHostName: String

  // /**
  //   * Builds the list of headers for the token request
  //   */
  // protected def buildTokenRequestHeaders(dockerHashContext: DockerHashContext): List[HttpHeader]
  
  // /* Methods that may be overridden by a subclass */

  // /**
  //   * service parameter that can be used in the token request
  //   * (e.g https://auth.docker.io/token?service=registry.docker.io&scope=repository:library/ubuntu:pull)
  //   */
  // protected def serviceName: Option[String] = None

  protected def regionId: String = "cn-shanghai"

  protected def endpointName: String = regionId

  protected def product: String = "cr"

  val regionPattern = """[^\s]+"""
  val hostPattern: Regex = s"""registry.($regionPattern).aliyuncs.com""".r

  // protected def domain: String = {
  //   s"$product.$regionId.aliyuncs.com"
  //}

  /**
    * Builds the token URI to be queried based on a DockerImageID
    */
  // protected def buildTokenRequestUri(dockerImageID: DockerImageIdentifierWithoutHash): String = {
  //   val service = serviceName map { name => s"service=$name&" } getOrElse ""
  //   s"https://$authorizationServerHostName/token?${service}scope=repository:${dockerImageID.nameWithDefaultRepository}:pull"
  // }

  // def getReginIdFromContext(context: DockerHashContext) = {
  //   val fullImageIdReg: Reg = s
  // }

  private def getManifest(/*notUsed: DockerHashResponse,*/ context: DockerHashContext) = {
    // val getCredentials = getAliyunProfileFromContext(context) match {

    // }



    val regionId = context.dockerImageID.host match {
      case Some(hostPattern(region)) => region
      case _ => "not match"
    }
    getAliyunCredentialFromContext(context) match {
      case None => manifestResponseHandler(new JavaHttpResponse, context)
      case Some(credentials: BasicCredentials) => {
        println("AliyunCr get manifest start")
        //DockerHashActor.logger.info(s"Attempting to pull {}", dockerCliKey.fullName)
        val endpoint = "cr." + regionId + ".aliyuncs.com"
        DefaultProfile.addEndpoint(regionId, AliyunCrAbstractFlow.ProductName, endpoint)
        
        // val profile: IClientProfile = getAliyunCredentialFromContext(context) match { 
        //   case Some(credentials: BasicCredentials) => DefaultProfile.getProfile("cn-shanghai", credentials.getAccessKeyId(), credentials.getAccessSecret())
        //   case Some(credentials: BasicSessionCredentials) => DefaultProfile.getProfile("cn-shanghai", credentials.getAccessKeyId(), credentials.getAccessSecret(), credentials.getSecurityToken())
        //   case None => 
        // }
        //val profile: IClientProfile = DefaultProfile.getProfile(defautProfle.getRegionId(), aliyunConfig.accessKeyId, aliyunConfig.accessKeySecret)
        val profile: IClientProfile = DefaultProfile.getProfile(regionId, credentials.getAccessKeyId(), credentials.getAccessKeySecret())
        val client: DefaultAcsClient = new DefaultAcsClient(profile)
        // 构造请求
        val request: GetRepoTagsRequest = new GetRepoTagsRequest()
        // 设置参数
        val dockerImageID = context.dockerImageID
        request.setRepoName(dockerImageID.image)
        dockerImageID.repository match {
          case Some(repository) => request.setRepoNamespace(repository)
          case _ => 
        }
        
        //request.setTag(dockerImageID.reference)
        // 发起请求
        val response = client.doAction(request)
        println(response)
        manifestResponseHandler(response, context)
      }
      case _ => manifestResponseHandler(new JavaHttpResponse, context)
    }
  }

  private def getAliyunCredentialFromContext(context: DockerHashContext): Option[Any] = {
    context.credentials find {_.isInstanceOf[BasicCredentials]/* || _.isInstanceOf[BasicCredentials]*/} match {
      case None => None
      case Some(cred: BasicCredentials) => Some(cred)
      //case Some(sCred: BasicSessionCredentials) => Some(sCred)
      case Some(_) => None
    }
  }

  // /**
  //   * Http method used for the manifest request
  //   */
  // // protected def manifestRequestHttpMethod: HttpMethod = HttpMethods.HEAD

  // /**
  //   * Generic method to build a flow that creates a request, sends it,
  //   * gets the result and performs some transformation on the result
  //   * @param requestBuilderFunction Function that builds a request from an input object of type A
  //   * @param responseHandlerFlow Function that transforms a response in an object of type Try[B]
  //   * @tparam A type of the input value of the flow
  //   * @tparam B type of the output value of the flow, if successful
  //   * @return flow
  //   */
  // private def requestTransformFlow[A, B](
  //                                         requestBuilderFunction: A => (HttpRequest, DockerHashContext),
  //                                         responseHandlerFlow: Flow[(HttpResponse, DockerHashContext), (Try[B], DockerHashContext), NotUsed]
  //                                       ) = GraphDSL.create() { implicit builder =>
  //   import GraphDSL.Implicits._

  //   // Build a request from an A
  //   val request = builder.add(Flow.fromFunction(requestBuilderFunction))

  //   // Perform http requests
  //   val http = builder.add(httpFlowWithRetry)
    
  //   // We got a HttpResponse, it might still be an unsuccessful response though.
  //   // Regardless send it through the responseHandlerFlow
  //   val finalResponse = http.out0.via(responseHandlerFlow)
    
  //   // The request could not be executed, map it to a failure message
  //   val failedRequest = http.out1.map({
  //     case (failure, commandContext) => throwableToFailureMessage(failure, commandContext)
  //   })

  //   // Partition processed responses
  //   val processResponse = builder.add(FlowUtils.fanOutTry[B, DockerHashContext])
    
  //   val successfulProcessing = processResponse.out0
  //   val failedProcessing = processResponse.out1.map(Function.tupled(throwableToFailureMessage))

  //   // Merge failures from failed response and failed parsing
  //   val mergeFailures = builder.add(Merge[(DockerHashResponse, DockerHashContext)](2))

  //   request ~> http.in
  //   finalResponse ~> processResponse.in

  //   failedRequest ~> mergeFailures.in(0)
  //   failedProcessing ~> mergeFailures.in(1)

  //   new FanOutShape2(request.in, successfulProcessing, mergeFailures.out)
  // }

  // /**
  //   * Builds the token request
  //   */
  // private def buildTokenRequest(dockerHashContext: DockerHashContext) = {
  //   val authorizationHeaders = buildTokenRequestHeaders(dockerHashContext)

  //   val tokenHttpRequest = HttpRequest(
  //     method = HttpMethods.GET,
  //     uri = buildTokenRequestUri(dockerHashContext.dockerImageID),
  //     headers = authorizationHeaders
  //   )

  //   (tokenHttpRequest, dockerHashContext)
  // }

  // /**
  //   * Parse the http response coming back from the token request and returns the body as JsObject
  //   */
  // private def tokenResponseHandler(response: HttpResponse, dockerHashContext: DockerHashContext) = response match {
  //   case httpResponse if httpResponse.status.isSuccess() =>
  //     Unmarshal(httpResponse.entity).to[String] map { body =>
  //       (
  //         Try(body.parseJson.asJsObject) flatMap extractToken,
  //         dockerHashContext
  //       )
  //     } recoverWith {
  //       case failure => Future.successful((Failure(failure), dockerHashContext))
  //     }
  //   case httpResponse =>
  //     Future.successful(
  //       (
  //         Failure(UnsuccessfulHttpResponseException(httpResponse)),
  //         dockerHashContext
  //       )
  //     )
  // }

  // /**
  //   * Extract the access token from the json body of the http response
  //   */
  // private def extractToken(jsObject: JsObject): Try[String] = {
  //   jsObject.fields.get("token") match {
  //     case Some(token: JsString) => Success(token.value)
  //     case Some(_) => Failure(new Exception("Token response contains a non-string token field"))
  //     case None => Failure(new Exception("Token response did not contain a token field"))
  //   }
  // }

  /**
    * Builds the manifest URI to be queried based on a DockerImageID
    */
  // private def buildManifestUri(dockerImageID: DockerImageIdentifierWithoutHash): String = {
  //   s"https://$registryHostName/v2/${dockerImageID.nameWithDefaultRepository}/manifests/${dockerImageID.reference}"
  // }

  /**
    * Builds the manifest http request
    */
  // private def buildManifestRequest(token: Option[String], dockerHashContext: DockerHashContext) = {
  //   val manifestRequest = token match {
  //     case Some(authToken) =>
  //       HttpRequest(
  //         method = manifestRequestHttpMethod,
  //         uri = buildManifestUri(dockerHashContext.dockerImageID),
  //         headers = scala.collection.immutable.Seq(AcceptHeader, Authorization(OAuth2BearerToken(authToken)))
  //       )
  //     case None =>
  //       HttpRequest(
  //         method = manifestRequestHttpMethod,
  //         uri = buildManifestUri(dockerHashContext.dockerImageID)
  //       )
  //   }

  //   (manifestRequest, dockerHashContext)
  // }

  private def matchTag(jsObject: JsObject, dockerHashContext: DockerHashContext): Boolean = {
    val tag = dockerHashContext.dockerImageID.reference
    jsObject.fields.get("tag") match {
      case Some(tagObj: JsString) if tagObj.value == tag => true
      case _ => false 
    }
  }

  private def extractDigestFromBody(jsObject: JsObject, dockerHashContext: DockerHashContext): Try[DockerHashResponse]= {
    val dataObj = jsObject.fields.get("data").get.asJsObject().convertTo[Map[String, JsValue]]
    val tags = dataObj.get("tags").get.convertTo[Seq[JsObject]]

    println(tags)
    println(tags.getClass())

    tags find {matchTag(_, dockerHashContext)} match {
      case Some(tagObj) =>
        println(tagObj.fields.get("digest"))
        tagObj.fields.get("digest") match {
          case Some(digest: JsString) => 
            DockerHashResult.fromString(AliyunCrAbstractFlow.HashAlg + ":" + digest.value) map { DockerHashSuccessResponse(_, dockerHashContext.request)}
          case Some(_) => Failure(new Exception("Manifest response contains a non-string digest field"))
          case None => Failure(new Exception("Manifest response did not contain a digest field"))
        }
      case None => 
        Failure(new Exception("Manifest response did not contain a digest field"))
    }
  }

  /**
    * Handles the http response from the manifest request
    */
  private def manifestResponseHandler(response: JavaHttpResponse, dockerHashContext: DockerHashContext) = response match {
    case httpResponse if httpResponse.isSuccess() =>
      Unmarshal(HttpEntity(httpResponse.getHttpContent)).to[String] map { body =>
        (
          Try(body.parseJson.asJsObject) flatMap {
            case jsObject => extractDigestFromBody(jsObject, dockerHashContext)
            },
          dockerHashContext
        )
      } recoverWith {
        case failure => Future.successful((Failure(failure), dockerHashContext))
      }
    case httpResponse =>
      Future.successful(
        (
          Failure(UnsuccessfulAliyunHttpResponseException(httpResponse)),
          dockerHashContext
        )
      )
    
    // case httpResponse if httpResponse.isSuccess() =>
    //   (
    //     extractContentDigest(response.headers, dockerHashContext),
    //     dockerHashContext
    //     )
    // case httpResponse => 
    //   (
    //     Failure(UnsuccessfulHttpResponseException(httpResponse)),
    //     dockerHashContext
    //  )
  }
  
  /**
    * Extracts the digest from the response headers
    */
  // private def extractContentDigest(headers: Seq[HttpHeader], dockerHashContext: DockerHashContext) = {
  //   headers find { _.is(DigestHeaderName) } match {
  //     case Some(digestHeader) =>
  //       DockerHashResult.fromString(digestHeader.value()) map { DockerHashSuccessResponse(_, dockerHashContext.request)}
  //     case None => Failure(new Exception("Cannot find digest header"))
  //   }
  // }

  /**
    * Maps a throwable into a DockerHashFailedResponse
    */
  private def throwableToFailureMessage(throwable: Throwable, dockerHashContext: DockerHashContext): (DockerHashResponse, DockerHashContext) = {
    throwable match {
      case UnsuccessfulHttpResponseException(httpResponse) => httpFailureMapper(httpResponse, dockerHashContext)
      case other => (DockerHashFailedResponse(other, dockerHashContext.request), dockerHashContext)
    }
  }

  /**
    * Maps a failed http response to a DockerHashFailureResponse
    */
  protected def httpFailureMapper(response: HttpResponse, dockerHashContext: DockerHashContext) = {
    response.status match {
      case StatusCodes.NotFound => (DockerHashNotFound(dockerHashContext.request), dockerHashContext)
      case StatusCodes.Unauthorized => (DockerHashUnauthorized(dockerHashContext.request), dockerHashContext)
      case other => (
        DockerHashFailedResponse(
          new Exception(s"Docker hash lookup failed with code ${other.intValue()}. ${other.defaultMessage()}"),
          dockerHashContext.request
        ),
        dockerHashContext
        )
    }
  }
}
