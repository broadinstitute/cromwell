package cromwell.services.instrumentation.impl.noop

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

object NoopInstrumentationServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(
    new NoopInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor)
  )
}

/**
  * Actor that ignores every InstrumentationServiceMessage - This is the default implementation of this service
  */
class NoopInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends Actor {
  override def receive = {
    case InstrumentationServiceMessage(_) => // Drop all messages
    case ShutdownCommand => context stop self
  }
}
