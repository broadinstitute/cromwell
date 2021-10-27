package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val azureKeyVaultName = configurationDescriptor.getOptBackendConfString("azure-keyvault-name")
  val azureB2CTokenSecretName = configurationDescriptor.getOptBackendConfString("azure-token-secret")
}
