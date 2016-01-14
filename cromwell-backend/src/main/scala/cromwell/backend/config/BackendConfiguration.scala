package cromwell.backend.config

import com.typesafe.config.ConfigFactory

import scala.collection.JavaConversions._

/**
  * Defines a backend configuration.
  * @param name Backend name.
  * @param initClass Initialization class for the specific backend.
  * @param `type` Backend type.
  */
case class BackendConfigurationEntry(name: String, initClass: String, `type`: String)

object BackendConfiguration {
  val config = ConfigFactory.load()
  val backendCfg = config.getConfig("backend")
  val defaultBackend = backendCfg.getString("default")
  val backendProviders = backendCfg.getConfigList("providers")
  val backendList = backendProviders.map(entry =>
    BackendConfigurationEntry(entry.getString("name"), entry.getString("initClass"), entry.getString("type"))).toList

  def apply(): BackendConfiguration = new BackendConfiguration(backendList, defaultBackend)
}

/**
  * Retrieves backend configuration.
  * @param backendList List of backend entries in configuration file.
  * @param defaultBackend Backend name to be used as default.
  */
class BackendConfiguration(backendList: List[BackendConfigurationEntry], defaultBackend: String) {

  /**
    * Gets default backend configuration.
    * @return Backend configuration.
    */
  def getDefaultBackend(): BackendConfigurationEntry = {
    backendList.filter(entry => entry.name.equals(defaultBackend)).head
  }

  /**
    * Gets all backend configurations from config file.
    * @return
    */
  def getAllBackendConfigurations(): List[BackendConfigurationEntry] = backendList

}
