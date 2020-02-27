package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor.DeleteMetadataAction
import cromwell.services.metadata.impl.MetadataDatabaseAccess

import scala.util.{Failure, Success}
import scala.concurrent.duration._

class DeleteMetadataActor(metadataDeletionConfig: MetadataDeletionConfig) extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  implicit val ec = context.dispatcher

  metadataDeletionConfig.intervalOpt.map { interval =>
    context.system.scheduler.schedule(5.minutes, interval, self, DeleteMetadataAction)
  }

  override def receive: Receive = {
    case DeleteMetadataAction =>
      val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(metadataDeletionConfig.delayAfterWorkflowCompletion.toSeconds)
      val workflowIdsForMetadataDeletionFuture = queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay,
        metadataDeletionConfig.batchSize
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
    case unexpected => log.warning(s"Programmer error: unexpected message $unexpected received from $sender")
  }
}

object DeleteMetadataActor {

  def props(metadataDeletionConfig: MetadataDeletionConfig) = Props(new DeleteMetadataActor(metadataDeletionConfig))

  case object DeleteMetadataAction

}
