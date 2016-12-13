package cromwell.services

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{Actor, ActorInitializationException, ActorLogging, ActorRef, OneForOneStrategy, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import com.typesafe.config.{Config, ConfigFactory, ConfigObject}
import net.ceedubs.ficus.Ficus._
import scala.collection.JavaConverters._

object ServiceRegistryActor {
  case class ServiceRegistryFailure(serviceName: String)

  trait ServiceRegistryMessage {
    def serviceName: String
  }

  def props(config: Config) = Props(new ServiceRegistryActor(serviceNameToPropsMap(config))).withDispatcher(ServiceDispatcher)

  // To enable testing, this lets us override a config value with a Props of our choice:
  def props(config: Config, overrides: Map[String, Props]) = {
    val fromConfig = serviceNameToPropsMap(config).filterNot { case (name: String, _: Props) => overrides.keys.toList.contains(name) }
    Props(new ServiceRegistryActor(fromConfig ++ overrides)).withDispatcher(ServiceDispatcher)
  }

  def serviceNameToPropsMap(globalConfig: Config): Map[String, Props] = {
    val serviceNamesToConfigStanzas = globalConfig.getObject("services").entrySet.asScala.map(x => x.getKey -> x.getValue).toMap
    serviceNamesToConfigStanzas map {
      case (serviceName, config: ConfigObject) => serviceName -> serviceProps(serviceName, globalConfig, config.toConfig)
      case (serviceName, _) => throw new Exception(s"Invalid configuration for service $serviceName")
    }
  }

  private def serviceProps(serviceName: String, globalConfig: Config, serviceStanza: Config): Props = {
    val serviceConfigStanza = serviceStanza.as[Option[Config]]("config").getOrElse(ConfigFactory.parseString(""))
    val className = serviceStanza.as[Option[String]]("class").getOrElse(
      throw new IllegalArgumentException(s"Invalid configuration for service $serviceName: missing 'class' definition")
    )
    try {
      Props.create(Class.forName(className), serviceConfigStanza, globalConfig)
    } catch {
      case e: ClassNotFoundException => throw new RuntimeException(
        s"Class $className for service $serviceName cannot be found in the class path.", e
      )
    }
  }
}

class ServiceRegistryActor(serviceProps: Map[String, Props]) extends Actor with ActorLogging {
  import ServiceRegistryActor._

  val services: Map[String, ActorRef] = serviceProps map {
    case (name, props) => name -> context.actorOf(props, name)
  }

  def receive = {
    case msg: ServiceRegistryMessage =>
      services.get(msg.serviceName) match {
        case Some(ref) => ref.tell(msg, sender)
        case None =>
          log.error("Received ServiceRegistryMessage requesting service '{}' for which no service is configured.  Message: {}", msg.serviceName, msg)
          sender ! ServiceRegistryFailure(msg.serviceName)
      }
    case fool =>
      log.error("Received message which is not a ServiceRegistryMessage: {}", fool)
      sender ! ServiceRegistryFailure("Message is not a ServiceRegistryMessage: " + fool)
  }

  /**
    * Set the supervision strategy such that any of the individual service actors fails to initialize that we'll pass
    * the error up the chain
    */
  override val supervisorStrategy = OneForOneStrategy() {
    case aie: ActorInitializationException => Escalate
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }
}

