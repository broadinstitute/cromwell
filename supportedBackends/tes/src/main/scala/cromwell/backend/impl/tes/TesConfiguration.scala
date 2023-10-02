package cromwell.backend.impl.tes

import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.retry.SimpleExponentialBackoff
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

class TesConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val endpointURL = configurationDescriptor.backendConfig.getString("endpoint")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val useBackendParameters =
    configurationDescriptor
      .backendConfig
      .as[Option[Boolean]](TesConfiguration.useBackendParametersKey)
      .getOrElse(false)

  val pollBackoff =
    configurationDescriptor
      .backendConfig
      .as[Option[Config]]("poll-backoff")
      .map(SimpleExponentialBackoff(_))
      .getOrElse(TesConfiguration.defaultPollBackoff)

  val executeOrRecoverBackoff =
    configurationDescriptor
      .backendConfig
      .as[Option[Config]]("execute-or-recover-backoff")
      .map(SimpleExponentialBackoff(_))
      .getOrElse(TesConfiguration.defaultExecOrRecoverBackoff)

  // Used for testing only. Include a bearer token for authenticating with the TES server
  final val bearerPrefix: String = "Bearer "
  val token: Option[String] = {
    configurationDescriptor.backendConfig.as[Option[String]]("bearer-token").map { t =>
      if (!t.startsWith(bearerPrefix))
        s"${bearerPrefix}${t}"
      else t
    }
  }
}

object TesConfiguration {
  final val useBackendParametersKey = "use_tes_11_preview_backend_parameters"

  final val defaultPollBackoff = SimpleExponentialBackoff(
    initialInterval = 10 seconds,
    maxInterval = 5 minutes,
    multiplier = 1.1
  )

  final val defaultExecOrRecoverBackoff = SimpleExponentialBackoff(
    initialInterval = 3 seconds,
    maxInterval = 30 seconds,
    multiplier = 1.1
  )
}
