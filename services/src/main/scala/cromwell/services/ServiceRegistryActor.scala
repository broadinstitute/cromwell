package cromwell.services

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.{Config, ConfigFactory, ConfigObject}
import lenthall.config.ScalaConfig._

import scala.collection.JavaConverters._

object ServiceRegistryActor {
  case class ServiceRegistryFailure(serviceName: String)
  trait ServiceRegistryMessage {
    def serviceName: String
  }

  def props(config: Config) = {
    Props(ServiceRegistryActor(config))
  }

  def serviceNameToPropsMap(globalConfig: Config): Map[String, Props] = {
    val serviceNamesToConfigStanzas = globalConfig.getObject("services").entrySet.asScala.map(x => x.getKey -> x.getValue).toMap
    serviceNamesToConfigStanzas map {
      case (serviceName, config: ConfigObject) =>
        System.out.println(s"Make me a service! $serviceName")
        serviceName -> serviceProps(serviceName, globalConfig, config.toConfig)
      case (serviceName, _) => throw new Exception(s"Invalid configuration for service $serviceName")
    }
  }

  private def serviceProps(serviceName: String, globalConfig: Config, serviceStanza: Config): Props = {
    val serviceConfigStanza = serviceStanza.getConfigOr("config", ConfigFactory.parseString(""))
    val className = serviceStanza.getStringOr(
      "class",
      throw new IllegalArgumentException(s"Invalid configuration for service $serviceName: missing 'class' definition")
    )
    Props.create(Class.forName(className), serviceConfigStanza, globalConfig)
  }
}

case class ServiceRegistryActor(globalConfig: Config) extends Actor with ActorLogging {
  import ServiceRegistryActor._

  val services: Map[String, ActorRef] = serviceNameToPropsMap(globalConfig) map {
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
}
