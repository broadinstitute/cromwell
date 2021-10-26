package cromwell.backend.impl.tes

import com.azure.identity.DefaultAzureCredentialBuilder
import com.azure.security.keyvault.secrets.{SecretClient, SecretClientBuilder}
import common.validation.ErrorOr
import common.validation.ErrorOr._

class AzureKeyVaultClient(client: SecretClient) {
  def getSecret(secretName: String): ErrorOr[String] = {
    ErrorOr(client.getSecret(secretName)).map(_.getValue)
  }
}

object AzureKeyVaultClient {
  def apply(vaultName: String, identityClientId: Option[String]): ErrorOr[AzureKeyVaultClient] = ErrorOr {
    val defaultCreds = identityClientId.map(identityId =>
      new DefaultAzureCredentialBuilder().managedIdentityClientId(identityId)
    ).getOrElse(
      new DefaultAzureCredentialBuilder()
    ).build()

    val client = new SecretClientBuilder()
      .vaultUrl(s"https://${vaultName}.vault.azure.net")
      .credential(defaultCreds)
      .buildClient()

    new AzureKeyVaultClient(client)
  }
}