package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val azureKeyVaultName = configurationDescriptor.backendConfig.getString("azure-keyvault-name")
  val azureB2CTokenSecretName = configurationDescriptor.backendConfig.getString("azure-token-secret")
}
