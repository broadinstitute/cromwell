package cromwell.services

import akka.actor.{Actor, ActorRef}
import com.typesafe.config.Config
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

class NooPServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor {
  override def receive = { case ShutdownCommand =>
    context stop self
  }
}
