package cromwell.util.docker

import cromwell.util.google.GoogleCredentialFactory

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
  * Connects to a Docker Registry using the API.
  */
trait DockerRegistryApiClient {

  /** The context used for Future.flatMap, etc. */
  implicit def executionContext: ExecutionContext

  /**
    * Returns a Docker hashable from a Docker identifier.
    * @param id The identifier.
    * @return The hashable docker value.
    */
  final def getDockerHashable(id: DockerIdentifier): Future[DockerHashable] = {
    id match {
      case DockerDigestIdentifier(_, digest, _) => Future.successful(DockerDigestHashable(digest))
      case tagId: DockerTagIdentifier => getManifest(tagId) fallbackTo getImageId(tagId)
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

object DockerLogin {
  lazy val GcrLogin = {
    // TODO: Are we supposed to be storing a refresh token / access token in a credential store?
    // TODO: This works for single use for now... will need to keep access token cached/refreshed, but not too rapidly?
    // TODO: http://stackoverflow.com/q/28272849/3320205

    val credential = GoogleCredentialFactory.fromAuthScheme
    // Our "user" auth scheme is authorizing, while the "service" auth scheme does not populate the credentials' tokens.
    if (credential.getAccessToken == null) credential.refreshToken()
    require(credential.getAccessToken != null, "GCR token unavailable, even after refresh")
    DockerLogin("_token", credential.getAccessToken)
  }

  // TODO: For local testing only right now. Longer term, we will NOT test GCR on Travis.
  lazy val GcrLoginOption: Option[DockerLogin] =
    try {
      Option(GcrLogin)
    } catch {
    case ex: Exception =>
      Console.err.println("Unable to get Google login information.")
      ex.printStackTrace(Console.err)
      None
  }
}

/** A login and password combo. */
case class DockerLogin(username: String, password: String)

/** Thrown/returned when an expected header is not found during an HTTP request. */
class DockerHeaderNotFoundException(name: String) extends Exception(s"Response did not contain header $name")

/** Thrown/returned when a docker manifest is not found. */
class DockerManifestNotFoundException(cause: Throwable)
  extends Exception(s"Did not find a manifest for the docker tag.", cause)
