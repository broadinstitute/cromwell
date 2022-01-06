package cromwell.backend.impl.tes

import cromwell.backend.BackendConfigurationDescriptor

import net.ceedubs.ficus.Ficus._

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val useBackendParameters =
    configurationDescriptor
      .backendConfig
      .as[Option[Boolean]](TesConfiguration.useBackendParamtersKey)
      .getOrElse(false)
}

object TesConfiguration {
  final val useBackendParamtersKey = "use_tes_11_preview_backend_parameters"
}
