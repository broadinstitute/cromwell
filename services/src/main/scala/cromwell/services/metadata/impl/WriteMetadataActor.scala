package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataService.{MetadataPutAcknowledgement, MetadataPutFailed, PutMetadataAction}

import scala.util.{Failure, Success}

object WriteMetadataActor {
  def props() = Props(new WriteMetadataActor()).withDispatcher(ServiceDispatcher)
}

class WriteMetadataActor extends Actor with ActorLogging with MetadataDatabaseAccess with SingletonServicesStore {

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
