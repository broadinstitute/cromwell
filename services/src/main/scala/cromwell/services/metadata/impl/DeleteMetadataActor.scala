package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.impl.DeleteMetadataActor.DeleteMetadataAction

import scala.util.{Failure, Success}
import scala.concurrent.duration._

class DeleteMetadataActor(timeToWaitAfterWorkflowFinish: Duration) extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  implicit val ec = context.dispatcher

  override def receive: Receive = {
    case DeleteMetadataAction =>
      val workflowIdsForMetadataDeletionFuture = queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        OffsetDateTime.now().minusSeconds(timeToWaitAfterWorkflowFinish.toSeconds)
      )
      workflowIdsForMetadataDeletionFuture onComplete {
        case Success(workflowIds) =>
          workflowIds foreach { workflowIdStr =>
            deleteNonLabelMetadataEntriesForWorkflowAndUpdateArchiveStatus(WorkflowId.fromString(workflowIdStr), MetadataArchiveStatus.toDatabaseValue(ArchivedAndPurged)) onComplete {
              case Success(_) => log.info(s"Successfully deleted metadata for workflow $workflowIdStr")
              case Failure(ex) => log.error(ex, s"Cannot delete metadata for workflow $workflowIdStr")
            }
          }
        case Failure(ex) =>
          log.error(ex, "Cannot delete metadata: unable to query list of workflow ids for metadata deletion from metadata summary table.")
      }
    case unknown@_ => throw new RuntimeException(s"Programmer's error: unknown message $unknown received from $sender")
  }
}

object DeleteMetadataActor {

  def props(timeToWaitAfterWorkflowFinish: Duration) = Props(new DeleteMetadataActor(timeToWaitAfterWorkflowFinish))

  case object DeleteMetadataAction

}
