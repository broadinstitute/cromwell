package cloud.nio.impl.drs

import cats.syntax.validated._
import com.azure.core.credential.TokenRequestContext
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.google.auth.oauth2.{AccessToken, GoogleCredentials, OAuth2Credentials}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.Ficus._

import java.time.{Duration => jDuration}
import java.util.concurrent.TimeUnit
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

/**
  * This trait allows us to abstract away different token attainment strategies
  * for different cloud environments.
  **/
sealed trait DrsCredentials {
  def getAccessToken: ErrorOr[String]
}

/**
  * Strategy for obtaining an access token from an existing OAuth credential. This class
  * is designed for use within the Cromwell engine.
  */
case class GoogleOauthDrsCredentials(credentials: OAuth2Credentials, acceptableTTL: Duration) extends DrsCredentials {
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
  private final val UserInfoEmailScope = "https://www.googleapis.com/auth/userinfo.email"
  private final val UserInfoProfileScope = "https://www.googleapis.com/auth/userinfo.profile"

  def getAccessToken: ErrorOr[String] = {
    Try {
      val scopedCredentials = GoogleCredentials.getApplicationDefault().createScoped(UserInfoEmailScope, UserInfoProfileScope)
      scopedCredentials.refreshAccessToken().getTokenValue
    } match {
      case Success(null) => "null token value attempting to refresh access token".invalidNel
      case Success(value) => value.validNel
      case Failure(e) => s"Failed to refresh access token: ${e.getMessage}".invalidNel
    }
  }
}


case class AzureDrsCredentials(identityClientId: Option[String]) extends DrsCredentials {

  final val tokenAcquisitionTimeout = new jDuration(30, TimeUnit.SECONDS)

  val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  val tokenScope = "https://management.azure.com/.default"

  def tokenRequestContext: TokenRequestContext = {
    val trc = new TokenRequestContext()
    trc.addScopes(tokenScope)
    trc
  }

  def defaultCredentialBuilder: DefaultAzureCredentialBuilder =
    new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)

  def getAccessToken: ErrorOr[String] = {
    val credentials = identityClientId.foldLeft(defaultCredentialBuilder) {
      (builder, clientId) => builder.managedIdentityClientId(clientId)
    }.build()

    Try(
      credentials
        .getToken(tokenRequestContext)
        .block(tokenAcquisitionTimeout)
    ) match {
      case Success(token) => token.getToken.validNel
      case Failure(error) =>
        Option(error.getMessage)  // it's possible the message is null
          .getOrElse(error.toString)
          .invalidNel
    }
  }
}
