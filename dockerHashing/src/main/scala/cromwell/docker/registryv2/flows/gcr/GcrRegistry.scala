package cromwell.docker.registryv2.flows.gcr

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import org.http4s.headers.Authorization
import org.http4s.{AuthScheme, Uri}

import scala.concurrent.duration._

class GcrRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  private val AccessTokenAcceptableTTL = 1.minute
  
  override protected def buildTokenRequestUri(dockerImageID: DockerImageIdentifier): Uri = {
    val uri = super.buildTokenRequestUri(dockerImageID)
    uri.withPath(s"/v2${uri.path}")
  }

  /**
    * Builds the list of headers for the token request
    */
   override def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext) = {
    dockerInfoContext.credentials collect {
      case credentials: OAuth2Credentials => Authorization(org.http4s.Credentials.Token(AuthScheme.Bearer, freshAccessToken(credentials)))
    }
  }
  
  private def freshAccessToken(credential: OAuth2Credentials) = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken) = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(AccessTokenAcceptableTTL)
    }
    
    Option(credential.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) => accessToken.getTokenValue
      case _ =>
        credential.refresh()
        credential.getAccessToken.getTokenValue
    }
  }

  override def accepts(dockerImageIdentifier: DockerImageIdentifier) = {
    dockerImageIdentifier.hostAsString.contains("gcr.io")
  }
}
