package cromwell.docker.registryv2.flows.google

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import org.http4s.headers.Authorization
import org.http4s.{AuthScheme, Uri}

import scala.concurrent.duration._

class GoogleRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  private val AccessTokenAcceptableTTL = 1.minute

  def googleRegion(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host.flatMap(_.split("/").headOption).getOrElse("")

  override def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = googleRegion(
    dockerImageIdentifier
  )
  override def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String = googleRegion(
    dockerImageIdentifier
  )

  override protected def buildTokenRequestUri(dockerImageID: DockerImageIdentifier): Uri = {
    val uri = super.buildTokenRequestUri(dockerImageID)
    uri.withPath(s"/v2${uri.path}")
  }

  /**
    * Builds the list of headers for the token request
    */
  def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Authorization] =
    dockerInfoContext.credentials collect { case credentials: OAuth2Credentials =>
      Authorization(org.http4s.Credentials.Token(AuthScheme.Bearer, freshAccessToken(credentials)))
    }

  private def freshAccessToken(credential: OAuth2Credentials) = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken) =
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(AccessTokenAcceptableTTL)

    Option(credential.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) => accessToken.getTokenValue
      case _ =>
        credential.refresh()
        credential.getAccessToken.getTokenValue
    }
  }

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean =
    // Supports both GCR (Google Container Registry) and GAR (Google Artifact Registry).
    dockerImageIdentifier.hostAsString.contains("gcr.io") || dockerImageIdentifier.hostAsString.contains(
      "-docker.pkg.dev"
    )
}
