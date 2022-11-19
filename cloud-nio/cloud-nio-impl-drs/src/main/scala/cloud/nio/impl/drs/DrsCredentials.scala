package cloud.nio.impl.drs

import cats.syntax.validated._
import com.azure.core.credential.TokenRequestContext
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
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
