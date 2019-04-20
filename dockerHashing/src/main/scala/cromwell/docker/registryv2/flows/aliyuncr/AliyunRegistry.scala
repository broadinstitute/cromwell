package cromwell.docker.registryv2.flows.aliyunregistry


import spray.json._
//import java.util.concurrent.TimeoutException
//import akka.http.scaladsl.unmarshalling.Unmarshal
import cats.effect.{ IO}
//import cromwell.docker.DockerInfoActor._
import cromwell.docker._
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.client.Client
import org.http4s.Header
import akka.stream._

import cromwell.docker.DockerInfoActor._
import cromwell.docker.{ DockerHashResult}
import spray.json._
import DefaultJsonProtocol._


import scala.concurrent.duration._
import scala.util.{Failure, Success}
import scala.util.matching.Regex
import com.aliyuncs.auth.{AlibabaCloudCredentials, BasicCredentials, BasicSessionCredentials}
import com.aliyuncs.DefaultAcsClient
import com.aliyuncs.cr.model.v20160607.GetRepoTagsRequest
import com.aliyuncs.http.{HttpResponse => JavaHttpResponse}
import com.aliyuncs.profile.DefaultProfile
import com.aliyuncs.profile.IClientProfile
import spray.json.{JsObject, JsString, JsValue}

/**
  * A docker flow using the CLI to return docker hashes.
  */
class AliyunRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  //  implicit val cs = IO.contextShift(ec)
  //  implicit val timer = IO.timer(ec)

  val ProductName = "cr"
  val HashAlg = "sha256"
  lazy val firstLookupTimeout = 5.seconds
  val supportAliunCrRegion = List("cn-qingdao", "cn-beijing", "cn-zhangjiakou", "cn-huhehaote",
    "cn-hangzhou", "cn-shanghai", "cn-shenzhen", "cn-hongkong",
    "ap-northeast-1", "ap-southeast-1", "ap-southeast-2", "ap-southeast-3",
    "ap-southeast-5", "ap-south-1", "us-east-1", "us-west-1", "me-east-1",
    "eu-central-1").mkString("|")
  val validAliyunCrHosts: Regex = ("""registry.""" + s"""(?:($supportAliunCrRegion))""" + """.aliyuncs.com""").r

  def isValidAliyunCrHost(host: Option[String]): Boolean =
    host match {
      case Some(h) => {
        h match {
          case validAliyunCrHosts(_) => true
          case _ => false
        }
      }
      case _ => false
    }

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = isValidAliyunCrHost(dockerImageIdentifier.host)

  override protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO.pure(None)
  }

  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = ""
  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String = ""
  // Not used for now, same reason as above
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] = List.empty

  override protected def getDockerResponse(token: Option[String], dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[DockerInfoSuccessResponse] = {
    //getManifest(dockerInfoContext)
    //IO.fromFuture(getManifest(dockerInfoContext))
    getManifest(dockerInfoContext) match {
        case succ: DockerInfoSuccessResponse => IO(succ)
        case fail: DockerInfoFailedResponse => throw fail.failure
        case other => throw new Exception(s"Get manifest failed, $other")
    }
  }

  val regionPattern = """[^\s]+"""
  val hostPattern: Regex = s"""registry.($regionPattern).aliyuncs.com""".r

  private def getManifest(context: DockerInfoContext) = {

    val regionId = context.dockerImageID.host match {
      case Some(hostPattern(region)) => region
      case _ => ""
    }

    val endpoint = ProductName + "." + regionId + ".aliyuncs.com"
    DefaultProfile.addEndpoint(regionId, ProductName, endpoint)

    val profile: IClientProfile = getAliyunCredentialFromContext(context) match {
      case Some(cred: BasicCredentials) => DefaultProfile.getProfile(regionId, cred.getAccessKeyId(), cred.getAccessKeySecret())
      case Some(sCred: BasicSessionCredentials) => DefaultProfile.getProfile(regionId, sCred.getAccessKeyId(), sCred.getAccessKeySecret(), sCred.getSessionToken())
      case _ => DefaultProfile.getProfile()
    }

    val client: DefaultAcsClient = new DefaultAcsClient(profile)
    val request: GetRepoTagsRequest = new GetRepoTagsRequest()
    val dockerImageID = context.dockerImageID
    request.setRepoName(dockerImageID.image)
    dockerImageID.repository foreach { repository => request.setRepoNamespace(repository) }

    manifestResponseHandler(client.doAction(request), context)
  }

  private def getAliyunCredentialFromContext(context: DockerInfoContext): Option[AlibabaCloudCredentials] = {
    context.credentials find {
      _.isInstanceOf[AlibabaCloudCredentials]
    } match {
      case Some(cred: BasicCredentials) => Some(cred)
      case Some(sCred: BasicSessionCredentials) => Some(sCred)
      case _ => None
    }
  }

  private def matchTag(jsObject: JsObject, dockerHashContext: DockerInfoContext): Boolean = {
    val tag = dockerHashContext.dockerImageID.reference
    jsObject.fields.get("tag") match {
      case Some(tagObj: JsString) if tagObj.value == tag => true
      case _ => false
    }
  }

  private def extractDigestFromBody(jsObject: JsObject, dockerHashContext: DockerInfoContext): DockerInfoResponse = {
    val dataObj = jsObject.fields.get("data").get.asJsObject().convertTo[Map[String, JsValue]]
    val tags = dataObj.get("tags").get.convertTo[Seq[JsObject]]

    tags find { matchTag(_, dockerHashContext)} match {
      case Some(tagObj) =>
        tagObj.fields.get("digest") match {
          case Some(digest: JsString) =>
            DockerHashResult.fromString(HashAlg + ":" + digest.value) match {
              case Success(r) => DockerInfoSuccessResponse(DockerInformation(r, None), dockerHashContext.request)
              case Failure(t) => DockerInfoFailedResponse(t, dockerHashContext.request)
            }
          case Some(_) => DockerInfoFailedResponse((new Exception("Manifest response contains a non-string digest field")), dockerHashContext.request)
          case None => DockerInfoFailedResponse((new Exception("Manifest response did not contain a digest field")), dockerHashContext.request)
        }
      case None => DockerInfoFailedResponse((new Exception("Manifest response did not contain a digest field")), dockerHashContext.request)
    }
  }

  private def manifestResponseHandler(response: JavaHttpResponse, dockerHashContext: DockerInfoContext) = response match {
    case httpResponse if httpResponse.isSuccess() =>
      extractDigestFromBody(new String(response.getHttpContent).parseJson.asJsObject(), dockerHashContext)
    case httpResponse => DockerInfoFailedResponse((new Exception(s"Get manifest not success $httpResponse")), dockerHashContext.request)
  }
}

