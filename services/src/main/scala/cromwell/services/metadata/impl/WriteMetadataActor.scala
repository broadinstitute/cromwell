package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.database.CromwellDatabase
import cromwell.services.metadata.MetadataService.{MetadataPutAcknowledgement, MetadataPutFailed, PutMetadataAction}

import scala.util.{Failure, Success}

object WriteMetadataActor {
  def props() = Props(new WriteMetadataActor())
}

class WriteMetadataActor extends Actor with ActorLogging with MetadataDatabaseAccess with CromwellDatabase {

  implicit val ec = context.dispatcher

  def receive = {
    case action@PutMetadataAction(events) =>
      val sndr = sender()
      addMetadataEvents(events) onComplete {
        case Success(_) => sndr ! MetadataPutAcknowledgement(action)
        case Failure(t) =>
          val msg = MetadataPutFailed(action, t)
          log.error(t, "Sending {} failure message {}", sndr, msg)
          sndr ! msg
      }
  }
}
