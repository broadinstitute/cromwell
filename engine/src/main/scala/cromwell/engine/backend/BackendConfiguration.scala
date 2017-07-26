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
    s"""|The configuration value 'backend.default' has been replaced with a list of values 'backend.allowed'.
        |The first value in the list is the default. Please update your configuration and restart.
        |""".stripMargin)

  lazy val AllowedBackendNames = backendConfig.getStringList("allowed").asScala.toList
  lazy val DefaultBackendName = AllowedBackendNames.headOption.getOrElse(
    throw new IllegalArgumentException("Specified allowed backends must not be an empty list."))
  private lazy val BackendProviders = backendConfig.getConfig("providers")

  lazy val AllBackendEntries: List[BackendConfigurationEntry] = AllowedBackendNames map { backendName =>
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

  /**
    * The global singleton `BackendConfiguration` for all main code. Some test specs creates their own instances
    * of the trait `BackendConfiguration` by indirectly invoking `fromConfig`.
    */
  val Global = fromConfig(ConfigFactory.load)
}
