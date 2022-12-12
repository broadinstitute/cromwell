package cromwell.backend.google.pipelines.batch

import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor

class GcpBatchConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val root: String = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig: Option[Config] = configurationDescriptor.backendRuntimeAttributesConfig

}
