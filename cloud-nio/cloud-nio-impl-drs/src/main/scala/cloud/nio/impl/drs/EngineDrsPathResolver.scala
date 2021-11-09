package cloud.nio.impl.drs

import com.google.auth.oauth2.{AccessToken, OAuth2Credentials}
import common.validation.ErrorOr.ErrorOr

import scala.concurrent.duration._
import cats.syntax.validated._
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.security.keyvault.secrets.SecretClientBuilder
import com.google.cloud.NoCredentials
import common.validation.ErrorOr


case class EngineDrsPathResolver(drsConfig: DrsConfig,
                                 accessTokenAcceptableTTL: Duration,
                                 authCredentials: OAuth2Credentials,
                                )
  extends DrsPathResolver(drsConfig, retryInternally = false) {

  //Based on method from GoogleRegistry
  override def getAccessToken: ErrorOr[String] = {
    def accessTokenTTLIsAcceptable(accessToken: AccessToken): Boolean = {
      (accessToken.getExpirationTime.getTime - System.currentTimeMillis()).millis.gteq(accessTokenAcceptableTTL)
    }

    authCredentials match {
      case _: NoCredentials => AzureKeyVaultClient.getToken("jdewar-0cb51953-vault", "b2c-token", None)
      case creds: OAuth2Credentials => Option(creds.getAccessToken) match {
        case Some(accessToken) if accessTokenTTLIsAcceptable(accessToken) =>
          accessToken.getTokenValue.validNel
        case _ =>
          creds.refresh()
          Option(creds.getAccessToken.getTokenValue) match {
            case Some(accessToken) => accessToken.validNel
            case None => "Could not refresh access token".invalidNel
          }
      }
    }
  }
}

object AzureKeyVaultClient {
  def getToken(vaultName: String, secretName: String, identityClientId: Option[String]): ErrorOr[String] = ErrorOr {
    val defaultCreds = identityClientId.map(identityId =>
      new DefaultAzureCredentialBuilder().managedIdentityClientId(identityId)
    ).getOrElse(
      new DefaultAzureCredentialBuilder()
    ).build()

    val client = new SecretClientBuilder()
      .vaultUrl(s"https://${vaultName}.vault.azure.net")
      .credential(defaultCreds)
      .buildClient()

    client.getSecret(secretName).getValue
  }
}