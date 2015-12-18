package cromwell.util.docker

import scala.concurrent.{ExecutionContext, Future}

/**
  * Connects to a Docker Registry using the API.
  */
trait DockerRegistryApiClient {

  /** The context used for Future.flatMap, etc. */
  implicit def executionContext: ExecutionContext

  /**
    * Parses the identifier string using the default identifier parser, then returns the docker hash for the identifier,
    * either from the hash, or by contacting the identifier's docker registry.
    *
    * @param id The identifier string.
    * @return The docker hash.
    */
  final def getDockerHash(id: String): Future[DockerHash] = getDockerHash(DockerIdentifierParser.Default.parse(id))

  /**
    * Returns a docker hash for the identifier, either from the hash, or by contacting the identifier's docker registry.
    *
    * @param id The identifier.
    * @return The docker hash.
    */
  final def getDockerHash(id: DockerIdentifier): Future[DockerHash] = {
    for {
      hashable <- getDockerHashable(id)
      hash <- Future.fromTry(hashable.dockerHash)
    } yield hash
  }

  /**
    * Returns a Docker hashable from a Docker identifier.
    * @param id The identifier.
    * @return The hashable docker value.
    */
  final def getDockerHashable(id: DockerIdentifier): Future[DockerHashable] = {
    id match {
      case DockerDigestIdentifier(_, digest, _) => Future.successful(DockerDigestHashable(digest))
      // http://stackoverflow.com/questions/25946942/what-is-the-use-of-scala-future-fallbackto
      case tagId: DockerTagIdentifier => getManifest(tagId) recoverWith { case _ => getImageId(tagId) }
    }
  }

  /**
    * Returns a docker image id using the V1 API.
    */
  final def getImageId(tagId: DockerTagIdentifier): Future[DockerImageId] = {
    for {
      token <- getV1Token(tagId)
      imageId <- getImageId(tagId, token)
    } yield imageId
  }

  /**
    * Returns a docker manifest using the V2 API.
    */
  final def getManifest(tagId: DockerTagIdentifier): Future[DockerManifest] = {
    val futureManifestResponse = for {
      tokenRequest <- getV2TokenRequest(tagId)
      tokenResponse <- getV2TokenResponse(tagId, tokenRequest)
      manifestResponse <- getManifest(tagId, tokenResponse.token)
    } yield manifestResponse

    futureManifestResponse
  }

  /** Makes an HTTP request to the docker registry V1 API to retrieve the docker V1 token. */
  protected def getV1Token(tagId: DockerTagIdentifier): Future[DockerV1Token]

  /** Makes an HTTP request to the docker registry V1 API to retrieve a docker image ID. */
  protected def getImageId(tagId: DockerTagIdentifier, dockerV1Token: DockerV1Token): Future[DockerImageId]

  /** Makes an HTTP request to the docker registry V2 API to retrieve the request information for a V2 token. */
  protected def getV2TokenRequest(tagId: DockerTagIdentifier): Future[DockerV2TokenRequest]

  /** Makes an HTTP request to the docker registry V2 API to retrieve a docker registry V2 token. */
  protected def getV2TokenResponse(tagId: DockerTagIdentifier,
                                   tokenRequest: DockerV2TokenRequest): Future[DockerV2TokenResponse]

  /** Makes an HTTP request to the docker registry V2 API to retrieve a docker manifest. */
  protected def getManifest(tagId: DockerTagIdentifier, dockerV2Token: DockerV2Token): Future[DockerManifest]
}

/** Thrown/returned when an expected header is not found during an HTTP request. */
class DockerHeaderNotFoundException(name: String) extends Exception(s"Response did not contain header $name")

/** Thrown/returned when a docker manifest is not found. */
class DockerManifestNotFoundException(cause: Throwable)
  extends Exception(s"Did not find a manifest for the docker tag", cause)
