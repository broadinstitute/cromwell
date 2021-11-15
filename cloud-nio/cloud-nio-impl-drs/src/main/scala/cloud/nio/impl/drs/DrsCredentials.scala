package cloud.nio.impl.drs

import cats.syntax.validated._
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.security.keyvault.secrets.{SecretClient, SecretClientBuilder}
import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import com.typesafe.config.Config
import common.validation.ErrorOr
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

/**
  * This trait allows us to abstract away different token attainment strategies
  * for different cloud environments.
  **/
sealed trait DrsCredentials {
  def getAccessToken: ErrorOr[String]
}

case class GoogleDrsCredentials(credentials: OAuth2Credentials, acceptableTTL: Duration) extends DrsCredentials {
  //Based on method from GoogleRegistry
  def getAccessToken: ErrorOr[String] = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(acceptableTTL)
    }

    Option(credentials.getAccessToken) match {
      case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) =>
        accessToken.getTokenValue.validNel
      case _ =>
        credentials.refresh()
        Option(credentials.getAccessToken.getTokenValue) match {
          case Some(accessToken) => accessToken.validNel
          case None => "Could not refresh access token".invalidNel
        }
    }
  }
}

object GoogleDrsCredentials {
  def apply(credentials: OAuth2Credentials, config: Config): GoogleDrsCredentials =
    GoogleDrsCredentials(credentials, config.as[FiniteDuration]("access-token-acceptable-ttl"))
}

case class AzureDrsCredentials(identityClientId: Option[String], vaultName: String, secretName: String) extends DrsCredentials {

  lazy val secretClient: ErrorOr[SecretClient] = ErrorOr {
    val defaultCreds = identityClientId.map(identityId =>
      new DefaultAzureCredentialBuilder().managedIdentityClientId(identityId)
    ).getOrElse(
      new DefaultAzureCredentialBuilder()
    ).build()

    new SecretClientBuilder()
      .vaultUrl(s"https://${vaultName}.vault.azure.net")
      .credential(defaultCreds)
      .buildClient()
  }

  def getAccessToken: ErrorOr[String] = secretClient.map(_.getSecret(secretName).getValue)
}
