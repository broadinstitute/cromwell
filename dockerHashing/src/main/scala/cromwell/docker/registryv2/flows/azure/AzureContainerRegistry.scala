package cromwell.docker.registryv2.flows.azure

import cats.data.Validated.{Invalid, Valid}
import cats.effect.IO
import com.typesafe.scalalogging.LazyLogging
import common.validation.ErrorOr.ErrorOr
import cromwell.cloudsupport.azure.AzureCredentials
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.{Header, Request, Response, Status}
import cromwell.docker.registryv2.flows.azure.AzureContainerRegistry.domain
import org.http4s.circe.jsonOf
import org.http4s.client.Client
import io.circe.generic.auto._
import org.http4s._

class AzureContainerRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) with LazyLogging {

  /**
    * (e.g registry-1.docker.io)
    */
  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host.getOrElse("")

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean =
    dockerImageIdentifier.hostAsString.contains(domain)

  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host.getOrElse("")

  /**
    * In Azure, service name does not exist at the registry level, it varies per repo, e.g. `terrabatchdev.azurecr.io`
    */
  override def serviceName: Option[String] =
    throw new Exception("ACR service name is host of user-defined registry, must derive from `DockerImageIdentifier`")

  /**
    * Builds the list of headers for the token request
    */
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] =
    List(contentTypeHeader)

  private val contentTypeHeader: Header = {
    import org.http4s.headers.`Content-Type`
    import org.http4s.MediaType

    `Content-Type`(MediaType.application.`x-www-form-urlencoded`)
  }

  private def getRefreshToken(authServerHostname: String, defaultAccessToken: String): IO[Request[IO]] = {
    import org.http4s.Uri.{Authority, Scheme}
    import org.http4s.client.dsl.io._
    import org.http4s._

    val uri = Uri.apply(
      scheme = Option(Scheme.https),
      authority = Option(Authority(host = Uri.RegName(authServerHostname))),
      path = "/oauth2/exchange",
      query = Query.empty
    )

    org.http4s.Method.POST(
      UrlForm(
        "service" -> authServerHostname,
        "access_token" -> defaultAccessToken,
        "grant_type" -> "access_token"
      ),
      uri,
      List(contentTypeHeader): _*
    )
  }

  /*
  Unlike other repositories, Azure reserves `GET /oauth2/token` for Basic Authentication [0]
  In order to use Oauth we must `POST /oauth2/token` [1]

  [0] https://github.com/Azure/acr/blob/main/docs/Token-BasicAuth.md#using-the-token-api
  [1] https://github.com/Azure/acr/blob/main/docs/AAD-OAuth.md#calling-post-oauth2token-to-get-an-acr-access-token
   */
  private def getDockerAccessToken(hostname: String, repository: String, refreshToken: String): IO[Request[IO]] = {
    import org.http4s.Uri.{Authority, Scheme}
    import org.http4s.client.dsl.io._
    import org.http4s._

    val uri = Uri.apply(
      scheme = Option(Scheme.https),
      authority = Option(Authority(host = Uri.RegName(hostname))),
      path = "/oauth2/token",
      query = Query.empty
    )

    org.http4s.Method.POST(
      UrlForm(
        // Tricky behavior - invalid `repository` values return a 200 with a valid-looking token.
        // However, the token will cause 401s on all subsequent requests.
        "scope" -> s"repository:$repository:pull",
        "service" -> hostname,
        "refresh_token" -> refreshToken,
        "grant_type" -> "refresh_token"
      ),
      uri,
      List(contentTypeHeader): _*
    )
  }

  override protected def getToken(
    dockerInfoContext: DockerInfoContext
  )(implicit client: Client[IO]): IO[Option[String]] = {
    val hostname = authorizationServerHostName(dockerInfoContext.dockerImageID)
    val maybeAadAccessToken: ErrorOr[String] =
      AzureCredentials.getAccessToken(None) // AAD token suitable for get-refresh-token request
    val repository = dockerInfoContext.dockerImageID.image // ACR uses what we think of image name, as the repository

    // Top-level flow: AAD access token -> refresh token -> ACR access token
    maybeAadAccessToken match {
      case Valid(accessToken) =>
        (for {
          refreshToken <- executeRequest(getRefreshToken(hostname, accessToken), parseRefreshToken)
          dockerToken <- executeRequest(getDockerAccessToken(hostname, repository, refreshToken), parseAccessToken)
        } yield dockerToken).map(Option.apply)
      case Invalid(errors) =>
        IO.raiseError(
          new Exception(s"Could not obtain AAD token to exchange for ACR refresh token. Error(s): ${errors}")
        )
    }
  }

  implicit val refreshTokenDecoder: EntityDecoder[IO, AcrRefreshToken] = jsonOf[IO, AcrRefreshToken]
  implicit val accessTokenDecoder: EntityDecoder[IO, AcrAccessToken] = jsonOf[IO, AcrAccessToken]

  private def parseRefreshToken(response: Response[IO]): IO[String] = response match {
    case Status.Successful(r) => r.as[AcrRefreshToken].map(_.refresh_token)
    case r =>
      r.as[String]
        .flatMap(b => IO.raiseError(new Exception(s"Request failed with status ${r.status.code} and body $b")))
  }

  private def parseAccessToken(response: Response[IO]): IO[String] = response match {
    case Status.Successful(r) => r.as[AcrAccessToken].map(_.access_token)
    case r =>
      r.as[String]
        .flatMap(b => IO.raiseError(new Exception(s"Request failed with status ${r.status.code} and body $b")))
  }

}

object AzureContainerRegistry {

  def domain: String = "azurecr.io"

}
