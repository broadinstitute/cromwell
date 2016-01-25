package cromwell.backend.config

import com.typesafe.config.ConfigFactory

import scala.collection.JavaConversions._

/**
  * Defines a backend configuration.
  * @param name Backend name.
  * @param initClass Initialization class for the specific backend.
  */
case class BackendConfigurationEntry(name: String, initClass: String)

object BackendConfiguration {
  val config = ConfigFactory.load()
  val backendCfg = config.getConfig("backend")
  val defaultBackend = backendCfg.getString("default")
  val backendProviders = backendCfg.getConfigList("providers")
  val backendList = backendProviders.map(entry =>
    BackendConfigurationEntry(entry.getString("name"), entry.getString("initClass"))).toList

  def apply(): BackendConfiguration = new BackendConfiguration(backendList, defaultBackend)
}

/**
  * Retrieves backend configuration.
  * @param backendList List of backend entries in configuration file.
  * @param defaultBackend Backend name to be used as default.
  */
class BackendConfiguration(backendList: List[BackendConfigurationEntry], defaultBackend: String) {
  /**
    * Gets default backend configuration. There will be always just one default backend defined in configuration file.
    * It lookup for the backend definition which contains the name defined in 'default' entry in backend configuration.
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
