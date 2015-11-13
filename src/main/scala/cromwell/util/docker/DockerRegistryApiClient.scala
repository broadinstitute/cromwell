package cromwell.util.docker

import scala.concurrent.{ExecutionContext, Future}

trait DockerRegistryApiClient {

  implicit def executionContext: ExecutionContext

  final def getDockerHashable(id: DockerIdentifier, login: Option[DockerLogin]): Future[DockerHashable] = {
    id match {
      case DockerDigestIdentifier(_, digest, _) => Future.successful(DockerDigestHashable(digest))
      case tagId: DockerTagIdentifier => getManifest(tagId, login) fallbackTo getImageId(tagId, login)
    }
  }

  final def getImageId(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerImageId] = {
    for {
      token <- getV1Token(tagId, login)
      imageId <- getImageId(tagId, token)
    } yield imageId
  }

  final def getManifest(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerManifest] = {
    val futureManifestResponse = for {
      tokenRequest <- getV2TokenRequest(tagId)
      tokenResponse <- getV2TokenResponse(tagId, tokenRequest, login)
      manifestResponse <- getManifest(tagId, tokenResponse.token)
    } yield manifestResponse

    futureManifestResponse
  }

  protected def getV1Token(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerV1Token]

  protected def getImageId(tagId: DockerTagIdentifier, dockerV1Token: DockerV1Token): Future[DockerImageId]

  protected def getV2TokenRequest(tagId: DockerTagIdentifier): Future[DockerV2TokenRequest]

  protected def getV2TokenResponse(tagId: DockerTagIdentifier, tokenRequest: DockerV2TokenRequest,
                                   login: Option[DockerLogin]): Future[DockerV2TokenResponse]

  protected def getManifest(tagId: DockerTagIdentifier, dockerV2Token: DockerV2Token): Future[DockerManifest]
}

case class DockerLogin(username: String, password: String)

class DockerHeaderNotFoundException(name: String) extends Exception(s"Response did not contain header $name")

class DockerManifestNotFoundException(cause: Throwable)
  extends Exception(s"Did not find a manifest for the docker tag.", cause)
