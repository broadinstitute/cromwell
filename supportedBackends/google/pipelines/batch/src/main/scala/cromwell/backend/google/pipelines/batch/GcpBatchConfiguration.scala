package cromwell.backend.google.pipelines.batch

import cromwell.backend.BackendConfigurationDescriptor

class GcpBatchConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val root = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig

}
