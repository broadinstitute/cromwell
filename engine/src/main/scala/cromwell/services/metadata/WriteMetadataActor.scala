package cromwell.services.metadata

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.engine.db.DataAccess
import cromwell.services.MetadataServiceActor.{MetadataPutAcknowledgement, MetadataPutFailed, PutMetadataAction}

import scala.util.{Failure, Success}

object WriteMetadataActor {
  def props() = Props(new WriteMetadataActor())
}

class WriteMetadataActor extends Actor with ActorLogging {

  val dataAccess = DataAccess.globalDataAccess

  implicit val ec = context.dispatcher

  def receive = {
    case action@PutMetadataAction(events) =>
      val sndr = sender()
      dataAccess.addMetadataEvents(events) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
  }
}
