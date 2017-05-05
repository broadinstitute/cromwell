package cromwell.engine.backend

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.{BackendConfigurationDescriptor, BackendLifecycleActorFactory}
import net.ceedubs.ficus.Ficus._
import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

case class BackendConfigurationEntry(name: String, lifecycleActorFactoryClass: String, config: Config) {
  def asBackendLifecycleActorFactory: Try[BackendLifecycleActorFactory] = Try {
    Class.forName(lifecycleActorFactoryClass)
         .getConstructor(classOf[String], classOf[BackendConfigurationDescriptor])
         .newInstance(name, asBackendConfigurationDescriptor)
         .asInstanceOf[BackendLifecycleActorFactory]
  }

  def asBackendConfigurationDescriptor = BackendConfigurationDescriptor(config, ConfigFactory.load)
}

trait BackendConfiguration {
  protected def backendConfig: Config

  require(!backendConfig.hasPath("default"),
    s"""|The configuration value 'backend.default' has been replaced with a list of values 'backend.enabled'.
        |The first value in the list is the default. Please update your configuration and restart.
        |""".stripMargin)

  lazy val EnabledBackendNames = backendConfig.getStringList("enabled").asScala.toList
  lazy val DefaultBackendName = EnabledBackendNames.headOption.getOrElse(
    throw new IllegalArgumentException("Specified enabled backends must not be an empty list."))
  private lazy val BackendProviders = backendConfig.getConfig("providers")

  lazy val AllBackendEntries: List[BackendConfigurationEntry] = EnabledBackendNames map { backendName =>
    val entry = BackendProviders.getConfig(backendName)
    BackendConfigurationEntry(
      backendName,
      entry.getString("actor-factory"),
      entry.as[Option[Config]]("config").getOrElse(ConfigFactory.empty("empty"))
    )
  }

  def backendConfigurationDescriptor(backendName: String): Try[BackendConfigurationDescriptor] = {
    AllBackendEntries.collect({case entry if entry.name.equalsIgnoreCase(backendName) => entry.asBackendConfigurationDescriptor}).headOption match {
      case Some(descriptor) => Success(descriptor)
      case None => Failure(new Exception(s"invalid backend: $backendName"))
    }
  }
}

object BackendConfiguration {
  def fromConfig(outerConfig: Config): BackendConfiguration = {
    new BackendConfiguration {
      override protected def backendConfig: Config = outerConfig.getConfig("backend")
    }
  }

  val Global = fromConfig(ConfigFactory.load)
}
