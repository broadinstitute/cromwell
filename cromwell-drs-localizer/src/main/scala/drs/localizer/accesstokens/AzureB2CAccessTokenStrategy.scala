package drs.localizer.accesstokens

import cats.syntax.validated._
import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.security.keyvault.secrets.{SecretClient, SecretClientBuilder}
import common.validation.ErrorOr
import common.validation.ErrorOr.{ErrorOr,ShortCircuitingFlatMap}
import drs.localizer.CommandLineArguments

case class AzureB2CAccessTokenStrategy(commandLineArguments: CommandLineArguments) extends AccessTokenStrategy {
  override def getAccessToken(): ErrorOr[String] = {
    commandLineArguments match {
      case CommandLineArguments(_, _, _, _, Some(vault), Some(secret), clientId) =>
        AzureKeyVaultClient(vault, clientId) flatMap { _.getSecret(secret) }
      case invalid => s"Invalid command line arguments: $invalid".invalidNel
    }
  }
}

// Note: The AzureKeyVaultClient code here is basically a copy/paste of the code temporarily living in the TES backend.
// All the current KeyVault interaction in Cromwell is temporary, but the code in the TES backend might be even more
// temporary than this.
class AzureKeyVaultClient(client: SecretClient) {
  def getSecret(secretName: String): ErrorOr[String] = ErrorOr(client.getSecret(secretName).getValue)
}

object AzureKeyVaultClient {
  def apply(vaultName: String, identityClientId: Option[String]): ErrorOr[AzureKeyVaultClient] = ErrorOr {
    val credentialBuilder = identityClientId.foldLeft(new DefaultAzureCredentialBuilder()) {
      (builder, clientId) => builder.managedIdentityClientId(clientId)
    }

    val client = new SecretClientBuilder()
      .vaultUrl(s"https://${vaultName}.vault.azure.net")
      .credential(credentialBuilder.build())
      .buildClient()

    new AzureKeyVaultClient(client)
  }
}
