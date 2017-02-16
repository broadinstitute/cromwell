package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor


class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val root = configurationDescriptor.backendConfig.getString("root")
  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
}
