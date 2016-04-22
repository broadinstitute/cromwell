package cromwell.services

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.{ConfigFactory, Config, ConfigObject}
import cromwell.engine.db.DataAccess
import cromwell.services.ServiceRegistryActor._
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
      case (serviceName, config: ConfigObject) => serviceName -> serviceProps(serviceName, globalConfig, config.toConfig)
      case (serviceName, _) => throw new Exception(s"Invalid configuration for service $serviceName")
    }
  }

  private def serviceProps(serviceName: String, globalConfig: Config, serviceStanza: Config): Props = {
    val serviceConfigStanza = serviceStanza.getConfigOr("config", ConfigFactory.parseString(""))
    val className = serviceStanza.getStringOr(
      "class",
      throw new Exception(s"Invalid configuration for service $serviceName: missing 'class' definition")
    )
    Props.create(Class.forName(className), serviceConfigStanza, globalConfig)
  }
}

case class ServiceRegistryActor(globalConfig: Config) extends Actor {
  import ServiceRegistryActor._

  val services: Map[String, ActorRef] = serviceNameToPropsMap(globalConfig) map {
    case (name, props) => name -> context.actorOf(props, "ServiceRegistryActor")
  }

  def receive = {
    case msg: ServiceRegistryMessage =>
      services.get(msg.serviceName) match {
        case Some(ref) => ref.tell(msg, sender)
        case None => sender ! ServiceRegistryFailure(msg.serviceName)
      }
  }
}
