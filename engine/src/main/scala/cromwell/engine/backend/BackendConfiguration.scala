package cromwell.engine.backend

import akka.actor.ActorSystem
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.{BackendConfigurationDescriptor, BackendLifecycleActorFactory}

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import lenthall.config.ScalaConfig._

case class BackendConfigurationEntry(name: String, lifecycleActorFactoryClass: String, config: Config) {
  def asBackendLifecycleActorFactory(actorSystem: ActorSystem): BackendLifecycleActorFactory = {
    Class.forName(lifecycleActorFactoryClass)
         .getConstructor(classOf[BackendConfigurationDescriptor], classOf[ActorSystem])
         .newInstance(asBackendConfigurationDescriptor, actorSystem)
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
    BackendConfigurationEntry(
      backendName,
      entry.getString("actor-factory"),
      entry.getConfigOr("config")
    )
  }

  val DefaultBackendEntry: BackendConfigurationEntry = AllBackendEntries.find(_.name == DefaultBackendName) getOrElse {
    throw new IllegalArgumentException(s"Could not find specified default backend name '$DefaultBackendName' " +
      s"in '${BackendNames.mkString("', '")}'.")
  }

  def backendConfigurationDescriptor(backendName: String): Try[BackendConfigurationDescriptor] = {
    AllBackendEntries.collect({case entry if entry.name.equalsIgnoreCase(backendName) => entry.asBackendConfigurationDescriptor}).headOption match {
      case Some(descriptor) => Success(descriptor)
      case None => Failure(new Exception(s"invalid backend: $backendName"))
    }
  }
}
