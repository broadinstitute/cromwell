package cromwell.services

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{Actor, ActorInitializationException, ActorLogging, ActorRef, OneForOneStrategy, Props}
import akka.routing.Listen
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory, ConfigObject}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._

object ServiceRegistryActor {
  case object NoopMessage
  case class ServiceRegistryFailure(serviceName: String)

  trait ServiceRegistryMessage {
    def serviceName: String
  }

  trait ListenToMessage extends ServiceRegistryMessage

  sealed trait ServiceRegistryMetaRequest
  case object RequestIoActorRef extends ServiceRegistryMetaRequest
  final case class IoActorRef(ioActor: ActorRef) extends ServiceRegistryMetaRequest
  case object NoIoActorRefAvailable

  def props(config: Config) = Props(new ServiceRegistryActor(config)).withDispatcher(ServiceDispatcher)

  // To enable testing, this lets us override a config value with a Props of our choice:
  def props(config: Config, overrides: Map[String, Props]) = {
    Props(new ServiceRegistryActor(config) {
      override def serviceProps = super.serviceProps ++ overrides
    }).withDispatcher(ServiceDispatcher)
  }

  def serviceNameToPropsMap(globalConfig: Config, registryActor: ActorRef): Map[String, Props] = {
    val serviceNamesToConfigStanzas = globalConfig.getObject("services").entrySet.asScala.map(x => x.getKey -> x.getValue).toMap
    serviceNamesToConfigStanzas map {
      case (serviceName, config: ConfigObject) => serviceName -> serviceProps(serviceName, globalConfig, config.toConfig, registryActor)
      case (serviceName, _) => throw new Exception(s"Invalid configuration for service $serviceName")
    }
  }

  private def serviceProps(serviceName: String, globalConfig: Config, serviceStanza: Config, registryActor: ActorRef): Props = {
    val serviceConfigStanza = serviceStanza.as[Option[Config]]("config").getOrElse(ConfigFactory.parseString(""))

    val dispatcher = serviceStanza.as[Option[String]]("dispatcher").getOrElse(ServiceDispatcher)
    val className = serviceStanza.as[Option[String]]("class").getOrElse(
      throw new IllegalArgumentException(s"Invalid configuration for service $serviceName: missing 'class' definition")
    )

    try {
      Props.create(Class.forName(className), serviceConfigStanza, globalConfig, registryActor).withDispatcher(dispatcher)
    } catch {
      case e: ClassNotFoundException => throw new RuntimeException(
        s"Class $className for service $serviceName cannot be found in the class path.", e
      )
    }
  }
}

class ServiceRegistryActor(globalConfig: Config) extends Actor with ActorLogging with GracefulShutdownHelper {
  import ServiceRegistryActor._

  def serviceProps = serviceNameToPropsMap(globalConfig, self)

  // When the IO actor starts up, it can register itself here for other service registry actors to make use of it
  var ioActor: Option[ActorRef] = None

  val services: Map[String, ActorRef] = serviceProps map {
    case (name, props) => name -> context.actorOf(props, name)
  }
  
  private def transform(message: Any, from: ActorRef): Any = message match {
    case _: ListenToMessage => Listen(from)
    case _ => message
  }

  def receive = {
    case msg: ServiceRegistryMessage =>
      services.get(msg.serviceName) match {
        case Some(ref) => ref.tell(transform(msg, sender), sender)
        case None =>
          log.error("Received ServiceRegistryMessage requesting service '{}' for which no service is configured.  Message: {}", msg.serviceName, msg)
          sender ! ServiceRegistryFailure(msg.serviceName)
      }
    case meta: ServiceRegistryMetaRequest => meta match {
      case RequestIoActorRef => ioActor match {
        case Some(ref) => sender ! IoActorRef(ref)
        case None => sender ! NoIoActorRefAvailable
      }
      case IoActorRef(ref) =>
        if (ioActor.isEmpty) { ioActor = Option(ref) }
        else { log.error(s"Programmer Error: More than one IoActor is trying to register itself in the service registry ($ref will *NOT* replace the existing $ioActor)") }
    }
    case ShutdownCommand =>
      services.values.toList match {
        case Nil => context stop self
        case head :: tail => waitForActorsAndShutdown(NonEmptyList.of(head, tail: _*))
      }
    case NoopMessage => // Nothing to do - useful for streams that use this actor as a sink and want to send a message on completion
    case fool =>
      log.error("Received message which is not a ServiceRegistryMessage: {}", fool)
      sender ! ServiceRegistryFailure("Message is not a ServiceRegistryMessage: " + fool)
  }

  /**
    * Set the supervision strategy such that any of the individual service actors fails to initialize that we'll pass
    * the error up the chain
    */
  override val supervisorStrategy = OneForOneStrategy() {
    case _: ActorInitializationException => Escalate
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }
}
