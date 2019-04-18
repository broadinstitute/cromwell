package cromwell.docker.registryv2

import cats.data.NonEmptyList
import cats.effect.IO
import cats.syntax.either._
import common.validation.Validation._
import cromwell.docker.DockerInfoActor._
import cromwell.docker.registryv2.DockerRegistryV2Abstract._
import cromwell.docker.{DockerHashResult, DockerImageIdentifier, DockerRegistry, DockerRegistryConfig}
import io.circe.Decoder
import io.circe.generic.auto._
import org.http4s.Uri.{Authority, Scheme}
import org.http4s._
import org.http4s.circe._
import org.http4s.client.Client
import org.http4s.client.dsl.io._
import org.http4s.headers._
import org.http4s.util.CaseInsensitiveString

object DockerRegistryV2Abstract {
  implicit class EnhancedParseResult[A](val parseResult: ParseResult[A]) extends AnyVal {
    def unsafe(context: String): A =
      parseResult
        .leftMap(_.message)
        .leftMap(NonEmptyList.one)
        .unsafe(context)
  }

  val DigestHeaderName = CaseInsensitiveString("Docker-Content-Digest")
  val ManifestV2MediaType = "application/vnd.docker.distribution.manifest.v2+json"
  val ManifestListV2MediaType = "application/vnd.docker.distribution.manifest.list.v2+json"

  // If one of those fails it means someone changed one of the strings above to an invalid one.
  val ManifestV2MediaRange = MediaRange.parse(ManifestV2MediaType)
    .unsafe("Cannot parse invalid manifest v2 content type. Please report this error.")
  val ManifestListV2MediaRange = MediaRange.parse(ManifestListV2MediaType)
    .unsafe("Cannot parse invalid manifest list v2 content type. Please report this error.")
  val AcceptManifestV2Header = Accept.parse(ManifestV2MediaType)
    .unsafe("Cannot parse invalid manifest v2 Accept header. Please report this error.")

  implicit val entityManifestDecoder = jsonEntityDecoder[DockerManifest](ManifestV2MediaRange)
  implicit val entityManifestListDecoder = jsonEntityDecoder[DockerManifestList](ManifestListV2MediaRange)
  implicit val entityTokenDecoder = jsonOf[IO, DockerAccessToken]

  /**
    * Creates a Json decoder for the given type but allows a different media type than application/json
    * This is necessary because the docker registry API responds with an "application/vnd.docker.distribution.manifest.v2+json"
    * and not the traditional "application/json". Adapted from CirceInstances.jsonOf
    */
  private def jsonEntityDecoder[A](mediaRange: MediaRange)(implicit decoder: Decoder[A]): EntityDecoder[IO, A] = EntityDecoder.decodeBy[IO, A](mediaRange) { message =>
    CirceInstances.defaultJsonDecoder[IO]
      .decode(message, strict = false)
      .flatMap({ json =>
        decoder.decodeJson(json)
          .fold(
            failure =>
              DecodeResult.failure(
                InvalidMessageBodyFailure(s"Could not decode JSON: $json", Some(failure))),
            DecodeResult.success(_)
          )
      })
  }
  
  // Placeholder exceptions that can be carried through IO before being converted to a DockerInfoFailedResponse
  private class Unauthorized() extends Exception
  private class NotFound() extends Exception
}

/**
  * Implements logic to retrieve information (namely hash and compressed size) about a docker image
  * from a registry server conforming to the docker registry API v2 Specification:
  * https://docs.docker.com/registry/spec/api/
  */
abstract class DockerRegistryV2Abstract(override val config: DockerRegistryConfig) extends DockerRegistry {
  implicit val ec = config.executionContext
  implicit val cs = IO.contextShift(ec)
  implicit val timer = IO.timer(ec)

  /**
    * This is the main function. Given a docker context and an http client, retrieve information about the docker image.
    */
  def run(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[(DockerInfoResponse, DockerInfoContext)] = {
    val dockerResponse = for {
      token <- getToken(dockerInfoContext)
      dockerSuccessResponse <- getDockerResponse(token, dockerInfoContext)
    } yield dockerSuccessResponse

    // Always map failures to a DockerHashFailedResponse instead of letting the IO fail, this is important so that the stream
    // that is calling this function will not fail
    dockerResponse.attempt
      .map({
        case Left(_: Unauthorized) => DockerInfoUnauthorized(dockerInfoContext.request)
        case Left(_: NotFound) => DockerInfoNotFound(dockerInfoContext.request)
        case Left(failure) => DockerInfoFailedResponse(failure, dockerInfoContext.request)
        case Right(value) => value
      })
      .map(_ -> dockerInfoContext)
  }

  // Execute a request. No retries because they're expected to already be handled by the client
  private def executeRequest[A](request: IO[Request[IO]], handler: Response[IO] => IO[A])(implicit client: Client[IO]): IO[A] = {
    client.fetch[A](request)(handler)
  }

  /**
    * Obtain an authorization token for the subsequent manifest request. Return IO.pure(None) if no token is needed
    * @param dockerInfoContext context
    * @param client http client
    * @return an authorization token
    */
  protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    val request = buildTokenRequest(dockerInfoContext)
    executeRequest(request, tokenResponseHandler).map(Option.apply)
  }

  /**
    * Return a DockerInfoResponse. Makes appropriate requests using the auth token if provided
    * @param token auth token
    * @param dockerInfoContext context
    * @param client http client
    * @return docker info response
    */
  protected def getDockerResponse(token: Option[String], dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[DockerInfoSuccessResponse] = {
    val request = manifestRequest(token, dockerInfoContext.dockerImageID)
    executeRequest(request, handleManifestResponse(dockerInfoContext, token))
  }

  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  def accepts(dockerImageIdentifier: DockerImageIdentifier) = dockerImageIdentifier.host.contains(registryHostName(dockerImageIdentifier))

  /* Methods that must to be implemented by a subclass */


  private def hostnameFromDockerSpec(dockerImageIdentifier: DockerImageIdentifier) = dockerImageIdentifier.host.flatMap(_.split("/").headOption).getOrElse("")

  /**
    * By default the hostname is everything before the first slash except for registries like Docker Hub where the registry
    * and authorization server host names are not specified but well-known.
    */
  protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier) = hostnameFromDockerSpec(dockerImageIdentifier)

  protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier) = hostnameFromDockerSpec(dockerImageIdentifier)

  /**
    * Builds the list of headers for the token request
    */
  protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] = List.empty

  /* Methods that may be overridden by a subclass */

  /**
    * service parameter that can be used in the token request
    * (e.g https://auth.docker.io/token?service=registry.docker.io&scope=repository:library/ubuntu:pull)
    */
  protected def serviceName: Option[String] = None

  /**
    * Builds the token URI to be queried based on a DockerImageID
    */
  protected def buildTokenRequestUri(dockerImageID: DockerImageIdentifier): Uri = {
    val service = serviceName map { name => s"service=$name&" } getOrElse ""
    Uri.apply(
      scheme = Option(Scheme.https),
      authority = Option(Authority(host = Uri.RegName(authorizationServerHostName(dockerImageID)))),
      path = "/token",
      query = Query.fromString(s"${service}scope=repository:${dockerImageID.nameWithDefaultRepository}:pull")
    )
  }

  /**
    * Builds the token request
    */
  private def buildTokenRequest(dockerInfoContext: DockerInfoContext): IO[Request[IO]] = {
    val request = Method.GET(
      buildTokenRequestUri(dockerInfoContext.dockerImageID),
      buildTokenRequestHeaders(dockerInfoContext): _*
    )
    request
  }

  /**
    * Parse the http response coming back from the token request and returns the body as JsObject
    */
  private def tokenResponseHandler(response: Response[IO]): IO[String] = response match {
    case Status.Successful(r) => r.as[DockerAccessToken].map(_.token)
    case r => r.as[String]
      .flatMap(b => IO.raiseError(new Exception(s"Request failed with status ${r.status.code} and body $b")))
  }

  /**
    * Builds the manifest URI to be queried based on a DockerImageID
    */
  private def buildManifestUri(dockerImageID: DockerImageIdentifier): Uri = {
    Uri.apply(
      scheme = Option(Scheme.https),
      authority = Option(Authority(host = Uri.RegName(registryHostName(dockerImageID)))),
      path = s"/v2/${dockerImageID.nameWithDefaultRepository}/manifests/${dockerImageID.reference}"
    )
  }

  /**
    * Request to get the manifest, using the auth token if provided
    */
  private def manifestRequest(token: Option[String], imageId: DockerImageIdentifier): IO[Request[IO]] = {
    val authorizationHeader = token.map(t => Authorization(Credentials.Token(AuthScheme.Bearer, t)))
    val request = Method.GET(
      buildManifestUri(imageId),
      List(
        Option(AcceptManifestV2Header),
        authorizationHeader
      ).flatten: _*
    )
    request
  }

  private def handleManifestResponse(dockerInfoContext: DockerInfoContext, token: Option[String])(response: Response[IO])(implicit client: Client[IO]): IO[DockerInfoSuccessResponse] = {
    // Getting the manifest content is not strictly necessary but just a bonus to get the size. If it fails, log the error and return None
    def handleManifestAttempt(attempt: Either[Throwable, Option[DockerManifest]]): Option[DockerManifest] = attempt match {
      case Left(failure) =>
        logger.warn(s"Could not get manifest for ${dockerInfoContext.dockerImageID.fullName}", failure)
        None
      case Right(manifest) => manifest
    }

    for {
      hashResult <- getDigestFromResponse(response)
      maybeManifest <- parseManifest(dockerInfoContext.dockerImageID, token)(response).attempt.map(handleManifestAttempt)
    } yield DockerInfoSuccessResponse(DockerInformation(hashResult, maybeManifest.map(_.compressedSize).map(DockerSize.apply)), dockerInfoContext.request)
  }

  /**
    * Handles the http response from the manifest request
    * The response can be of 2 sorts:
    * - A manifest (https://docs.docker.com/registry/spec/manifest-v2-2/#image-manifest-field-descriptions)
    * - A manifest list which contains a list of pointers to other manifests (https://docs.docker.com/registry/spec/manifest-v2-2/#manifest-list)
    * 
    * When a manifest list is returned, we need to pick one of the manifest pointers and make another request for that manifest.
    * 
    * Because the different manifests in the list are (supposed to be) variations of the same image over different platforms,
    * we simply pick the first one here since we only care about the approximate size, and we don't expect it to change drastically
    * between platforms.
    * If that assumption turns out to be incorrect, a smarter decision may need to be made to choose the manifest to lookup. 
    */
  private def parseManifest(dockerImageIdentifier: DockerImageIdentifier, token: Option[String])(response: Response[IO])(implicit client: Client[IO]): IO[Option[DockerManifest]] = response match {
    case Status.Successful(r) if r.headers.exists(_.value.equalsIgnoreCase(ManifestV2MediaType)) =>
      r.as[DockerManifest].map(Option.apply)
    case Status.Successful(r) if r.headers.exists(_.value.equalsIgnoreCase(ManifestListV2MediaType)) =>
      r.as[DockerManifestList].flatMap({ dockerManifestList =>
        obtainManifestFromList(dockerManifestList, dockerImageIdentifier, token)
      })
    case _ => IO.pure(None)
  }

  private def obtainManifestFromList(dockerManifestList: DockerManifestList,
                                     dockerImageIdentifier: DockerImageIdentifier,
                                     token: Option[String])(implicit client: Client[IO]): IO[Option[DockerManifest]] = {
    dockerManifestList.manifests.headOption
      .map(_.digest)
      .map(dockerImageIdentifier.swapReference) match {
      case Some(identifierWithNewHash) =>
        val request = manifestRequest(token, identifierWithNewHash)
        executeRequest(request, parseManifest(dockerImageIdentifier, token))
      case None =>
        logger.error(s"The manifest list for ${dockerImageIdentifier.fullName} was empty. Cannot proceed to obtain the size of image")
        IO.pure(None)
    }
  }

  private def getDigestFromResponse(response: Response[IO]): IO[DockerHashResult] = response match {
    case Status.Successful(r) => extractDigestFromHeaders(r.headers)
    case Status.Unauthorized(_) => IO.raiseError(new Unauthorized)
    case Status.NotFound(_) => IO.raiseError(new NotFound)
    case failed => failed.as[String].flatMap(body => IO.raiseError(new Exception(s"Failed to get manifest: $body"))
    )
  }
  
  private def extractDigestFromHeaders(headers: Headers) = {
    headers.find(a => a.toRaw.name.equals(DigestHeaderName)) match {
      case Some(digest) => IO.fromEither(DockerHashResult.fromString(digest.value).toEither)
      case None => IO.raiseError(new Exception(s"Manifest response did not have a digest header"))
    }
  }
}
