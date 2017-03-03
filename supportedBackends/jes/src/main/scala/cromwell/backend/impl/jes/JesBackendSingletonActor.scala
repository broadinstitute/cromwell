package cromwell.backend.impl.jes

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.JesApiQueryManagerRequest
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive

final case class JesBackendSingletonActor(qps: Int Refined Positive) extends Actor with ActorLogging {

  val jesApiQueryManager = context.actorOf(JesApiQueryManager.props(qps))

  override def receive = {
    case apiQuery: JesApiQueryManagerRequest =>
      log.debug("Forwarding API query to JES API query manager actor")
      jesApiQueryManager.forward(apiQuery)
  }
}

object JesBackendSingletonActor {
  def props(qps: Int Refined Positive): Props = Props(JesBackendSingletonActor(qps)).withDispatcher(BackendDispatcher)
}
