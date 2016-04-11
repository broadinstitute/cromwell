package cromwell.engine.backend

import com.typesafe.config.{Config, ConfigFactory}

import scala.collection.JavaConverters._

/**
  * Defines a backend configuration.
  *
  * @param name      Backend name.
  * @param className Implementation class for the specific backend.
  */
case class BackendConfigurationEntry(name: String, className: String, config: Config)

object BackendConfiguration {
  private val BackendConfig = ConfigFactory.load().getConfig("backend")
  private val DefaultBackendName = BackendConfig.getString("default")
  private val BackendProvidersConfigs = BackendConfig.getConfigList("providers").asScala.toList

  val AllBackendEntries: List[BackendConfigurationEntry] = BackendProvidersConfigs map { entry =>
    BackendConfigurationEntry(entry.getString("name"), entry.getString("class"), entry.getConfig("config")) }

  val DefaultBackendEntry: BackendConfigurationEntry =
    AllBackendEntries find { _.name == DefaultBackendName } getOrElse { throw new IllegalArgumentException(s"Could not find specified default backend name '$DefaultBackendName'.") }
}
