package cromwell.docker.registryv2.flows.alibabacloudcrregistry

import cats.effect.IO
import com.aliyuncs.DefaultAcsClient
import com.aliyuncs.auth.{AlibabaCloudCredentials, BasicCredentials, BasicSessionCredentials}
import com.aliyuncs.cr.model.v20160607.GetRepoTagsRequest
import com.aliyuncs.profile.DefaultProfile
import com.aliyuncs.profile.IClientProfile
import cromwell.docker.DockerHashResult
import cromwell.docker.DockerInfoActor._
import cromwell.docker._
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.Header
import org.http4s.client.Client
import scala.util.matching.Regex
import scala.util.{Failure, Success, Try}
import spray.json.DefaultJsonProtocol._
import spray.json._

class AlibabaCloudCRRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  val ProductName = "cr"
  val HashAlg = "sha256"
  val regionPattern = """[^\s]+"""
  val validAlibabaCloudCRHosts: Regex = s"""registry.($regionPattern).aliyuncs.com""".r

  val validCrEndpoint: Regex = s"""cr.($regionPattern).aliyuncs.com""".r
  val validCrVpcEndpoint: Regex = s"""cr-vpc.($regionPattern).aliyuncs.com""".r


  def isValidAlibabaCloudCRHost(host: Option[String]): Boolean = {
    host.exists {
      _ match {
        case validAlibabaCloudCRHosts(_) => true
        case _ => false
      }
    }
  }

  def isValidAlibabaCloudCREndpoint(endpoint: String): Boolean = {
    endpoint match {
      case validCrEndpoint(_) => true
      case validCrVpcEndpoint(_) => true
      case _ => false
    }
  }


  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = isValidAlibabaCloudCRHost(dockerImageIdentifier.host)

  override protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO.pure(None)
  }

  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = ""
  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String = ""
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] = List.empty

  override protected def getDockerResponse(token: Option[String], dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[DockerInfoSuccessResponse] = {
    getManifest(dockerInfoContext) match {
        case success: DockerInfoSuccessResponse => IO(success)
        case fail: DockerInfoFailedResponse => IO.raiseError(new Exception(fail.reason))
        case other => IO.raiseError(new Exception(s"Get manifest failed, $other"))
    }
  }

  private def getManifest(context: DockerInfoContext): DockerInfoResponse = {

    val regionId = context.dockerImageID.host match {
      case Some(validAlibabaCloudCRHosts(region)) => region
      case _ => throw new Exception(s"The host ${context.dockerImageID.host} does not have the expected region id")
    }

    val defaultEndpoint = ProductName + "." + regionId + ".aliyuncs.com"
    val endpoint = getAliyunEndpointFromContext(context).getOrElse(defaultEndpoint)
    DefaultProfile.addEndpoint(regionId, ProductName, endpoint)

    val profile: IClientProfile = getAliyunCredentialFromContext(context) match {
      case Some(cred: BasicCredentials) => DefaultProfile.getProfile(regionId, cred.getAccessKeyId(), cred.getAccessKeySecret())
      case Some(sCred: BasicSessionCredentials) => DefaultProfile.getProfile(regionId, sCred.getAccessKeyId(), sCred.getAccessKeySecret(), sCred.getSessionToken())
      case _ => throw new Exception(s"Invalid credential from context, ${context}")
    }

    val client: DefaultAcsClient = new DefaultAcsClient(profile)
    val request: GetRepoTagsRequest = new GetRepoTagsRequest()
    val dockerImageID = context.dockerImageID
    request.setRepoName(dockerImageID.image)
    dockerImageID.repository foreach { repository => request.setRepoNamespace(repository) }

    manifestResponseHandler(client, request, context)
      .getOrElse(new Exception(s"handle response fail, please make sure the image id is correct: ${context.dockerImageID}")) match {
      case succ: DockerInfoSuccessResponse => succ
      case fail: DockerInfoFailedResponse => fail
      case ex: Exception => throw new Exception(s"Get AliyunCr manifest failed, ${ex.getMessage}")
    }
  }

  private[alibabacloudcrregistry] def getAliyunCredentialFromContext(context: DockerInfoContext): Option[AlibabaCloudCredentials] = {
    context.credentials find {
      _.isInstanceOf[AlibabaCloudCredentials]
    } match {
      case Some(cred: BasicCredentials) => Some(cred)
      case Some(sCred: BasicSessionCredentials) => Some(sCred)
      case _ => None
    }
  }

  //cr.cn-beijing.aliyuncs.com  or  cr-vpc.cn-beijing.aliyuncs.com
  private[alibabacloudcrregistry] def getAliyunEndpointFromContext(context: DockerInfoContext): Option[String] = {
    context.credentials collectFirst {
      case endpoint: String if (isValidAlibabaCloudCREndpoint(endpoint)) => endpoint
    }
  }

  private def matchTag(jsObject: JsObject, dockerHashContext: DockerInfoContext): Boolean = {
    val tag = dockerHashContext.dockerImageID.reference
    jsObject.fields.get("tag") match {
      case Some(tagObj: JsString) if tagObj.value == tag => true
      case _ => false
    }
  }

  private[alibabacloudcrregistry] def extractDigestFromBody(jsObject: JsObject, dockerHashContext: DockerInfoContext): DockerInfoResponse = {
    val tags = jsObject.fields.get("data") match {
      case Some(data) => data.asJsObject().convertTo[Map[String, JsValue]].get("tags") match {
        case Some(tag) => tag.convertTo[Seq[JsObject]]
        case None => throw new Exception(s"Manifest response did not contain a tags field, ${jsObject}")
      }
      case None => throw new Exception(s"Manifest response did not contain a data field, Please make sure the existence of image, ${jsObject}")
    }

    tags find { matchTag(_, dockerHashContext)} match {
      case Some(tagObj) =>
        tagObj.fields.get("digest") match {
          case Some(digest: JsString) =>
            DockerHashResult.fromString(HashAlg + ":" + digest.value) match {
              case Success(r) => DockerInfoSuccessResponse(DockerInformation(r, None), dockerHashContext.request)
              case Failure(t) => DockerInfoFailedResponse(t, dockerHashContext.request)
            }
          case Some(_) => DockerInfoFailedResponse((new Exception(s"Manifest response contains a non-string digest field, ${jsObject}")), dockerHashContext.request)
          case None => DockerInfoFailedResponse((new Exception(s"Manifest response did not contain a digest field, ${jsObject}")), dockerHashContext.request)
        }
      case None => DockerInfoFailedResponse((new Exception(s"Manifest response did not contain a expected tag: ${dockerHashContext.dockerImageID.reference}, ${jsObject}")), dockerHashContext.request)
    }
  }

  private def manifestResponseHandler(client: DefaultAcsClient, request: GetRepoTagsRequest, dockerHashContext: DockerInfoContext): Try[DockerInfoResponse] = {
    for {
      response <- Try(client.doAction(request))
      jsObj <- Try(if (response.isSuccess) response.getHttpContentString.parseJson.asJsObject()
                    else throw new Exception(s"Get manifest request not success: ${response}"))
      dockInfoRes <- Try(extractDigestFromBody(jsObj, dockerHashContext))
    } yield dockInfoRes
  }
}

