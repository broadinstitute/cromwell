package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService.{DeleteMetadataAction, DeleteMetadataFailedResponse, DeleteMetadataSuccessfulResponse}

import scala.util.{Failure, Success}

class DeleteMetadataActor extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  implicit val ec = context.dispatcher

  override def receive: Receive = {
    case action@DeleteMetadataAction(workflowId, replyTo, maxAttempts) =>
      val deleteFromDbFuture = deleteNonLabelMetadataEntriesForWorkflow(workflowId)
      deleteFromDbFuture onComplete {
        case Success(_) => replyTo ! DeleteMetadataSuccessfulResponse(workflowId)
        case Failure(ex) if maxAttempts > 1 =>
          val remainingAttempts = action.maxAttempts - 1
          log.error(ex, s"Cannot delete metadata. Remaining number of attempts: $maxAttempts.")
          self ! action.copy(maxAttempts = remainingAttempts)
        case Failure(ex) =>
          log.error(ex, s"Cannot delete metadata.")
          replyTo ! DeleteMetadataFailedResponse(workflowId, ex)
      }
  }
}

object DeleteMetadataActor {

  def props() = Props(new DeleteMetadataActor)

}
