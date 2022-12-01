package cromwell.filesystems.blob

import cats.implicits.catsSyntaxValidatedId
import com.azure.core.credential.TokenRequestContext
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import common.validation.ErrorOr.ErrorOr

import scala.concurrent.duration._
import scala.jdk.DurationConverters._

import scala.util.{Failure, Success, Try}

/**
  * Strategy for obtaining an access token in an environment with available Azure identity.
  * If you need to disambiguate among multiple active user-assigned managed identities, pass
  * in the client id of the identity that should be used.
  */
case object AzureCredentials {

  final val tokenAcquisitionTimeout = 5.seconds

  val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  val tokenScope = "https://management.azure.com/.default"

  private def tokenRequestContext: TokenRequestContext = {
    val trc = new TokenRequestContext()
    trc.addScopes(tokenScope)
    trc
  }

  private def defaultCredentialBuilder: DefaultAzureCredentialBuilder =
    new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)

  def getAccessToken(identityClientId: Option[String]): ErrorOr[String] = {
    val credentials = identityClientId.foldLeft(defaultCredentialBuilder) {
      (builder, clientId) => builder.managedIdentityClientId(clientId)
    }.build()

    Try(
      credentials
        .getToken(tokenRequestContext)
        .block(tokenAcquisitionTimeout.toJava)
    ) match {
      case Success(null) => "null token value attempting to obtain access token".invalidNel
      case Success(token) => token.getToken.validNel
      case Failure(error) => s"Failed to refresh access token: ${error.getMessage}".invalidNel
    }
  }
}
