package cromwell.util.docker

import akka.actor.ActorSystem
import akka.event.Logging
import spray.client.pipelining._
import spray.http.HttpHeaders.{Authorization, `WWW-Authenticate`}
import spray.http._
import spray.httpx.UnsuccessfulResponseException
import spray.httpx.unmarshalling._

import scala.concurrent.Future

class SprayDockerRegistryApiClient()(implicit system: ActorSystem) extends DockerRegistryApiClient {

  implicit val executionContext = system.dispatcher

  import SprayDockerRegistryApiMarshalling._

  private val log = Logging(system, getClass)

  private def logSendReceive: SendReceive = {
    logRequest(log) ~> sendReceive ~> logResponse(log)
  }

  private def addDockerLogin(login: Option[DockerLogin]): RequestTransformer = {
    login match {
      case Some(userPass) => addCredentials(BasicHttpCredentials(userPass.username, userPass.password))
      case None => identity[HttpRequest]
    }
  }

  override def getV1Token(taggedIdentifier: DockerTagIdentifier) = {
    val pipeline: HttpRequest => Future[DockerV1Token] =
      addDockerLogin(taggedIdentifier.registry.loginProvider.dockerLogin) ~>
        // https://github.com/docker/docker-registry/issues/517#issuecomment-86704783
        addHeader("X-Docker-Token", "true") ~>
        logSendReceive ~>
        checkResponseSuccess ~>
        toDockerV1Token

    pipeline(Get(taggedIdentifier.imagesUri))
  }

  private def checkResponseSuccess(response: HttpResponse): HttpResponse = {
    if (response.status.isSuccess) {
      response
    } else {
      throw new UnsuccessfulResponseException(response)
    }
  }

  private[docker] def toDockerV1Token(response: HttpResponse): DockerV1Token = {
    val tokenHeaders = response.headers collect {
      // Header matching is only lower case: https://github.com/spray/spray/issues/844
      case HttpHeader("x-docker-token", token) => token
    }

    tokenHeaders.headOption.getOrElse(throw new DockerHeaderNotFoundException("X-Docker-Token"))
  }

  override def getImageId(taggedIdentifier: DockerTagIdentifier, dockerToken: DockerV1Token) = {
    val pipeline: HttpRequest => Future[DockerImageId] =
      addHeader(Authorization(GenericHttpCredentials("Token", dockerToken))) ~>
        logSendReceive ~>
        unmarshal[DockerImageId]

    pipeline(Get(taggedIdentifier.tagsUri))
  }

  override def getV2TokenRequest(taggedIdentifier: DockerTagIdentifier) = {
    val pipeline: HttpRequest => Future[DockerV2TokenRequest] =
      logSendReceive ~>
        toDockerV2TokenRequest

    pipeline(Get(taggedIdentifier.manifestUri))
  }

  private[docker] def toDockerV2TokenRequest(response: HttpResponse): DockerV2TokenRequest = {
    val dockerV2TokenRequests = response.headers collect {
      case `WWW-Authenticate`(challenge :: _) =>
        DockerV2TokenRequest(challenge.realm, challenge.params("service"), challenge.params.get("scope"))
    }

    dockerV2TokenRequests.headOption.getOrElse(throw new DockerHeaderNotFoundException(`WWW-Authenticate`.name))
  }

  override def getV2TokenResponse(taggedIdentifier: DockerTagIdentifier,
                                  tokenRequest: DockerV2TokenRequest) = {
    val pipeline: HttpRequest => Future[DockerV2TokenResponse] =
      addDockerLogin(taggedIdentifier.registry.loginProvider.dockerLogin) ~>
        logSendReceive ~>
        unmarshal[DockerV2TokenResponse]

    val queryMap = taggedIdentifier.tokenRequestQueryMap(tokenRequest)
    val uri = Uri(tokenRequest.realm).withScheme("https").withQuery(queryMap)

    pipeline(Get(uri))
  }

  override def getManifest(taggedIdentifier: DockerTagIdentifier, dockerToken: DockerV2Token) = {
    val pipeline: HttpRequest => Future[DockerManifest] =
      addCredentials(OAuth2BearerToken(dockerToken)) ~>
        logSendReceive ~>
        unmarshal[DockerManifest]

    pipeline(Get(taggedIdentifier.manifestUri)) recoverWith mapManifestResponseExceptions
  }

  private def mapManifestResponseExceptions[T]: PartialFunction[Throwable, Future[T]] = {
    case e: UnsuccessfulResponseException if e.response.status == StatusCodes.NotFound =>
      Future.failed(new DockerManifestNotFoundException(e))
  }
}
