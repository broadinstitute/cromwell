package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor
import net.ceedubs.ficus.Ficus._

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val azureKeyVaultName = configurationDescriptor.backendConfig.as[Option[String]]("azure-keyvault-name")
  val azureB2CTokenSecretName = configurationDescriptor.backendConfig.as[Option[String]]("azure-token-secret")
}
