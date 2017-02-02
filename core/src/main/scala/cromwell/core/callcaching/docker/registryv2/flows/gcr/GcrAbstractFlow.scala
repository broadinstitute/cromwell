package cromwell.core.callcaching.docker.registryv2.flows.gcr

import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.stream.ActorMaterializer
import com.google.api.client.auth.oauth2.Credential
import cromwell.core.callcaching.docker.DockerHashActor.DockerHashContext
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow
import cromwell.core.callcaching.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

abstract class GcrAbstractFlow(httpClientFlow: HttpDockerFlow, host: String)(implicit ec: ExecutionContext, materializer: ActorMaterializer) extends DockerRegistryV2AbstractFlow(httpClientFlow)(ec, materializer) {
  
  private val AccessTokenAcceptableTTL = 1.minute.toSeconds
  
  override val registryHostName = host
  override val authorizationServerHostName = s"$host/v2"
  
  /**
    * Builds the list of headers for the token request
    */
   def buildTokenRequestHeaders(dockerHashContext: DockerHashContext) = {
    dockerHashContext.credentials collect {
      case credential: Credential => Authorization(OAuth2BearerToken(freshAccessToken(credential)))
    }
  }
  
  private def freshAccessToken(credential: Credential) = {
    if (credential.getExpiresInSeconds < AccessTokenAcceptableTTL) {
      credential.refreshToken()
    }
    credential.getAccessToken
  }
}
