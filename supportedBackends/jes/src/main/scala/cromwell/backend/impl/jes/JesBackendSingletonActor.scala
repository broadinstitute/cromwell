package cromwell.backend.impl.jes

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.DoPoll

final case class JesBackendSingletonActor(jesConfiguration: JesConfiguration) extends Actor with ActorLogging {

  val pollingActor = context.actorOf(JesApiQueryManager.props(jesConfiguration))

  override def receive = {
    case poll: DoPoll =>
      log.debug("Forwarding status poll to JES polling actor")
      pollingActor.forward(poll)
  }
}

object JesBackendSingletonActor {
  def props(jesConfiguration: JesConfiguration): Props = Props(JesBackendSingletonActor(jesConfiguration))
}
