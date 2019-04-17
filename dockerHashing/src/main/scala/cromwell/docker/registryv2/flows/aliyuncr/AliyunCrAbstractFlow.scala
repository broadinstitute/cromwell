package cromwell.docker.registryv2.flows.aliyuncr

import akka.NotUsed
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream._
import akka.stream.scaladsl.{Flow, GraphDSL, Merge}
import akka.stream.FlowShape
import cromwell.docker.DockerHashActor._
import cromwell.docker.registryv2.flows.aliyuncr.AliyunCrAbstractFlow._
import cromwell.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest
import cromwell.docker.registryv2.flows.FlowUtils
import cromwell.docker.{DockerFlow, DockerHashResult, DockerImageIdentifierWithoutHash}
import spray.json._
import DefaultJsonProtocol._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Try}
import scala.util.matching.Regex
import com.aliyuncs.auth.{AlibabaCloudCredentials, BasicCredentials, BasicSessionCredentials}
import com.aliyuncs.DefaultAcsClient
import com.aliyuncs.cr.model.v20160607.GetRepoTagsRequest
import com.aliyuncs.http.{HttpResponse => JavaHttpResponse}
import com.aliyuncs.profile.DefaultProfile
import com.aliyuncs.profile.IClientProfile

object AliyunCrAbstractFlow {
  type HttpDockerFlow = Flow[(HttpRequest, ContextWithRequest[DockerHashContext]), (Try[HttpResponse], ContextWithRequest[DockerHashContext]), NotUsed]

  val ProductName = "cr"

  val HashAlg = "sha256"

  case class UnsuccessfulHttpResponseException(httpResponse: HttpResponse) extends Exception
  case class UnsuccessfulAliyunHttpResponseException(httpResponse: JavaHttpResponse) extends Exception

}

/**
  * Implements logic to build a flow that can retrieve a docker hash 
  * from a Aliyun CR
  */
abstract class AliyunCrAbstractFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer) extends DockerFlow {
  
  /**
    * Builds a flow for this docker registry
    */
  def buildFlow() = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    val manifest = builder.add(
      Flow
        .fromFunction(getManifest _)
        .mapAsync(1)(identity)
    )

    // Partition processed responses
    val processResponse = builder.add(FlowUtils.fanOutTry[DockerHashResponse, DockerHashContext])

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

  val regionPattern = """[^\s]+"""
  val hostPattern: Regex = s"""registry.($regionPattern).aliyuncs.com""".r

  private def getManifest(context: DockerHashContext) = {

    val regionId = context.dockerImageID.host match {
      case Some(hostPattern(region)) => region
      case _ => ""
    }

    val endpoint = AliyunCrAbstractFlow.ProductName+ "." + regionId + ".aliyuncs.com"
    DefaultProfile.addEndpoint(regionId, AliyunCrAbstractFlow.ProductName, endpoint)

    val profile: IClientProfile = getAliyunCredentialFromContext(context) match {
      case Some(cred: BasicCredentials) => DefaultProfile.getProfile(regionId, cred.getAccessKeyId(), cred.getAccessKeySecret())
      case Some(sCred: BasicSessionCredentials) => DefaultProfile.getProfile(regionId, sCred.getAccessKeyId(), sCred.getAccessKeySecret(), sCred.getSessionToken())
      case _ => DefaultProfile.getProfile()
    }

    val client: DefaultAcsClient = new DefaultAcsClient(profile)
    val request: GetRepoTagsRequest = new GetRepoTagsRequest()
    val dockerImageID = context.dockerImageID
    request.setRepoName(dockerImageID.image)
    dockerImageID.repository foreach {repository => request.setRepoNamespace(repository)}

    manifestResponseHandler(client.doAction(request), context)
  }

  private def getAliyunCredentialFromContext(context: DockerHashContext): Option[AlibabaCloudCredentials] = {
    context.credentials find {_.isInstanceOf[AlibabaCloudCredentials]} match {
      case Some(cred: BasicCredentials) => Some(cred)
      case Some(sCred: BasicSessionCredentials) => Some(sCred)
      case _ => None
    }
  }

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

    tags find {matchTag(_, dockerHashContext)} match {
      case Some(tagObj) =>
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

  }

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
