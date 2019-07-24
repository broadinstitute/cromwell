package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.impl.MetadataArchivalActor.CheckForArchival

/**
  * For use by the standard metadata service implementation. If enabled, will look for metadata events
  * with the key Archiveable and delete all rows corresponding to that workflow
  */
// FIXME: Instrumentation
class MetadataArchivalActor extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  // FIXME: determine scheme to set up polling interval

  override def preStart() = {
    self ! CheckForArchival
  }

  override def receive = {
    case CheckForArchival =>
      // FIXME: Check for archiveable tasks
      //  FIXME: For each task, archive
      // FIXME: when done, schedule a message to self
  }
}

object MetadataArchivalActor {
  final case object CheckForArchival
}