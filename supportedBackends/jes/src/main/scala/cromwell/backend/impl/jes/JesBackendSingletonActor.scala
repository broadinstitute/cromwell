package cromwell.backend.impl.jes

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.DoPoll

final case class JesBackendSingletonActor(qps: Int) extends Actor with ActorLogging {

  val pollingActor = context.actorOf(JesApiQueryManager.props(qps))

  override def receive = {
    case poll: DoPoll =>
      log.debug("Forwarding status poll to JES polling actor")
      pollingActor.forward(poll)
  }
}

object JesBackendSingletonActor {
  def props(qps: Int): Props = Props(JesBackendSingletonActor(qps)).withDispatcher(BackendDispatcher)
}
