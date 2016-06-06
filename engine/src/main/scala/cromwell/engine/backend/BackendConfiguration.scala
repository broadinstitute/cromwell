package cromwell.engine.backend

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.{BackendConfigurationDescriptor, BackendLifecycleActorFactory}

import scala.collection.JavaConverters._
import scala.util.{Try, Success, Failure}

case class BackendConfigurationEntry(name: String, lifecycleActorFactoryClass: String, config: Config) {
  def asBackendLifecycleActorFactory: BackendLifecycleActorFactory = {
    Class.forName(lifecycleActorFactoryClass)
         .getConstructor(classOf[Config])
         .newInstance(config)
         .asInstanceOf[BackendLifecycleActorFactory]
  }

  def asBackendConfigurationDescriptor = BackendConfigurationDescriptor(config, ConfigFactory.load)
}

object BackendConfiguration {
  private val BackendConfig = ConfigFactory.load.getConfig("backend")
  private val DefaultBackendName = BackendConfig.getString("default")
  private val BackendProviders = BackendConfig.getConfig("providers")
  private val BackendNames: Set[String] = BackendProviders.entrySet().asScala.map(_.getKey.split("\\.").toSeq.head).toSet

  val AllBackendEntries: List[BackendConfigurationEntry] = BackendNames.toList map { backendName =>
    val entry = BackendProviders.getConfig(backendName)
    // TODO PBE: "class" can be removed when switching to the shadow world; "actor-factory" renamed?
    BackendConfigurationEntry(
      backendName,
      entry.getString("actor-factory"),
      entry.getConfig("config")
    )
  }

  val DefaultBackendEntry: BackendConfigurationEntry = AllBackendEntries.find(_.name == DefaultBackendName) getOrElse {
    throw new IllegalArgumentException(s"Could not find specified default backend name '$DefaultBackendName'.")
  }

  def backendConfigurationDescriptor(backendName: String): Try[BackendConfigurationDescriptor] = {
    AllBackendEntries.collect({case entry if entry.name.equalsIgnoreCase(backendName) => entry.asBackendConfigurationDescriptor}).headOption match {
      case Some(descriptor) => Success(descriptor)
      case None => Failure(new Exception(s"invalid backend: $backendName"))
    }
  }
}
