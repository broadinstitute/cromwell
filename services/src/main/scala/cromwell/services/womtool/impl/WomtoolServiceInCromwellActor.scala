package cromwell.services.womtool.impl

import akka.actor.ActorRef
import com.typesafe.config.Config
import cromwell.services.womtool.WomtoolServiceActor
import cromwell.services.womtool.WomtoolServiceMessages.DescribeResult

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends WomtoolServiceActor {
  override def receive: Receive = {
    case _ =>
     sender ! DescribeResult(Math.random > 0.5, List("hello", "world"))
  }
}
