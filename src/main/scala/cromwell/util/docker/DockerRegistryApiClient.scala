package cromwell.util.docker

import scala.concurrent.{ExecutionContext, Future}

/**
  * Connects to a Docker Registry using the API.
  */
trait DockerRegistryApiClient {

  /** The context used for Future.flatMap, etc. */
  implicit def executionContext: ExecutionContext

  /**
    * Returns a Docker hashable from a Docker identifier.
    * @param id The identifier.
    * @param login The optional login information.
    * @return The hashable docker value.
    */
  final def getDockerHashable(id: DockerIdentifier, login: Option[DockerLogin]): Future[DockerHashable] = {
    id match {
      case DockerDigestIdentifier(_, digest, _) => Future.successful(DockerDigestHashable(digest))
      case tagId: DockerTagIdentifier => getManifest(tagId, login) fallbackTo getImageId(tagId, login)
    }
  }

  /**
    * Returns a docker image id using the V1 API.
    */
  final def getImageId(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerImageId] = {
    for {
      token <- getV1Token(tagId, login)
      imageId <- getImageId(tagId, token)
    } yield imageId
  }

  /**
    * Returns a docker manifest using the V2 API.
    */
  final def getManifest(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerManifest] = {
    val futureManifestResponse = for {
      tokenRequest <- getV2TokenRequest(tagId)
      tokenResponse <- getV2TokenResponse(tagId, tokenRequest, login)
      manifestResponse <- getManifest(tagId, tokenResponse.token)
    } yield manifestResponse

    futureManifestResponse
  }

  /** Makes an HTTP request to the docker registry V1 API to retrieve the docker V1 token. */
  protected def getV1Token(tagId: DockerTagIdentifier, login: Option[DockerLogin]): Future[DockerV1Token]

  /** Makes an HTTP request to the docker registry V1 API to retrieve a docker image ID. */
  protected def getImageId(tagId: DockerTagIdentifier, dockerV1Token: DockerV1Token): Future[DockerImageId]

  /** Makes an HTTP request to the docker registry V2 API to retrieve the request information for a V2 token. */
  protected def getV2TokenRequest(tagId: DockerTagIdentifier): Future[DockerV2TokenRequest]

  /** Makes an HTTP request to the docker registry V2 API to retrieve a docker registry V2 token. */
  protected def getV2TokenResponse(tagId: DockerTagIdentifier, tokenRequest: DockerV2TokenRequest,
                                   login: Option[DockerLogin]): Future[DockerV2TokenResponse]

  /** Makes an HTTP request to the docker registry V2 API to retrieve a docker manifest. */
  protected def getManifest(tagId: DockerTagIdentifier, dockerV2Token: DockerV2Token): Future[DockerManifest]
}

/** A login and password combo. */
case class DockerLogin(username: String, password: String)

/** Thrown/returned when an expected header is not found during an HTTP request. */
class DockerHeaderNotFoundException(name: String) extends Exception(s"Response did not contain header $name")

/** Thrown/returned when a docker manifest is not found. */
class DockerManifestNotFoundException(cause: Throwable)
  extends Exception(s"Did not find a manifest for the docker tag.", cause)
