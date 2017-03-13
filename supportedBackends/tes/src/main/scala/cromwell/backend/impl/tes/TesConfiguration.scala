package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeAttrsConfig = configurationDescriptor.backendConfig.getConfig("default-runtime-attributes")
}
