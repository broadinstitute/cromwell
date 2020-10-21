package cromwell.services.instrumentation.impl.selectivetsv

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.services.instrumentation.{CromwellBucket, CromwellIncrement}
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

object SelectiveTsvInstrumentationServiceActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new SelectiveTsvInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor))
}

/**
  * Actor that ignores every InstrumentationServiceMessage - This is the default implementation of this service
  */
class SelectiveTsvInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor {

  var statesMap: Map[String, Int] = Map.empty

  override def receive = {
    case InstrumentationServiceMessage(CromwellIncrement(CromwellBucket(_, path))) if path.last == "starting" =>
      val pathString = path.init.mkString(".")
      val current = statesMap.getOrElse(pathString, 0)
      statesMap = statesMap + (pathString -> current + 1)
      println(s"Incrementing: $path")
    case ShutdownCommand => context stop self
  }
}
