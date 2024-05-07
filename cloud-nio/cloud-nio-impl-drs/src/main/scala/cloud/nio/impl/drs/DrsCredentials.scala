package cloud.nio.impl.drs

import cats.syntax.validated._
import com.google.auth.oauth2.{AccessToken, GoogleCredentials, OAuth2Credentials}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.Ficus._
import cromwell.cloudsupport.azure.AzureCredentials

import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

/**
  * This trait allows us to abstract away different token attainment strategies
  * for different cloud environments.
  **/
trait DrsCredentials {
  def getAccessToken: ErrorOr[String]
}

/**
  * Strategy for obtaining an access token from an existing OAuth credential. This class
  * is designed for use within the Cromwell engine.
  */
case class GoogleOauthDrsCredentials(credentials: OAuth2Credentials, acceptableTTL: Duration) extends DrsCredentials {
  // Based on method from GoogleRegistry
  def getAccessToken: ErrorOr[String] = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean =
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(acceptableTTL)

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

object GoogleOauthDrsCredentials {
  def apply(credentials: OAuth2Credentials, config: Config): GoogleOauthDrsCredentials =
    GoogleOauthDrsCredentials(credentials, config.as[FiniteDuration]("access-token-acceptable-ttl"))
}

/**
  * Strategy for obtaining an access token from Google Application Default credentials that are assumed to already exist
  * in the environment. This class is designed for use by standalone executables running in environments
  * that have direct access to a Google identity (ex. CromwellDrsLocalizer).
  */
case object GoogleAppDefaultTokenStrategy extends DrsCredentials {
  final private val UserInfoEmailScope = "https://www.googleapis.com/auth/userinfo.email"
  final private val UserInfoProfileScope = "https://www.googleapis.com/auth/userinfo.profile"

  def getAccessToken: ErrorOr[String] =
    Try {
      val scopedCredentials =
        GoogleCredentials.getApplicationDefault().createScoped(UserInfoEmailScope, UserInfoProfileScope)
      scopedCredentials.refreshAccessToken().getTokenValue
    } match {
      case Success(null) => "null token value attempting to refresh access token".invalidNel
      case Success(value) => value.validNel
      case Failure(e) => s"Failed to refresh access token: ${e.getMessage}".invalidNel
    }
}

/**
  * Strategy for obtaining an access token in an environment with available Azure identity.
  * If you need to disambiguate among multiple active user-assigned managed identities, pass
  * in the client id of the identity that should be used.
  */
case class AzureDrsCredentials(identityClientId: Option[String] = None) extends DrsCredentials {
  def getAccessToken: ErrorOr[String] = AzureCredentials.getAccessToken(identityClientId)
}
