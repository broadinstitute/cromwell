package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor._
import cromwell.services.metadata.hybridcarbonite.NumberOfWorkflowsToDeleteMetadataMetricActor.{CalculateNumberOfWorkflowsToDeleteMetadataMetricValue, NumberOfWorkflowsToDeleteMetadataMetricValue}
import cromwell.services.metadata.impl.MetadataDatabaseAccess

import scala.util.{Failure, Success}

class DeleteMetadataActor(metadataDeletionConfig: MetadataDeletionConfig, override val serviceRegistryActor: ActorRef) extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  var activeConfig: ActiveMetadataDeletionConfig = _
  metadataDeletionConfig match {
    case a: ActiveMetadataDeletionConfig =>
      logger.info(s"Archived metadata deletion is configured to begin polling after ${a.initialDelay}, and then delete up to ${a.batchSize} workflows worth of metadata every ${a.interval} (for workflows which completed at least ${a.delayAfterWorkflowCompletion} ago).")
      activeConfig = a
      context.system.scheduler.schedule(a.initialDelay, a.interval, self, DeleteMetadataAction)
    case InactiveMetadataDeletionConfig =>
      logger.info("Archived metadata deletion is configured to be inactive.")
      context.stop(self)
  }

  val numOfWorkflowsToDeleteMetadataMetricActor = context.actorOf(NumberOfWorkflowsToDeleteMetadataMetricActor.props(serviceRegistryActor))

  override def receive: Receive = {
    case DeleteMetadataAction =>
      val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(activeConfig.delayAfterWorkflowCompletion.toSeconds)
      val workflowIdsForMetadataDeletionFuture = queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay,
        activeConfig.batchSize
      )
      workflowIdsForMetadataDeletionFuture onComplete {
        case Success(workflowIds) =>
          if (workflowIds.length < metadataDeletionConfig.batchSize) {
            numOfWorkflowsToDeleteMetadataMetricActor ! NumberOfWorkflowsToDeleteMetadataMetricValue(workflowIds.length.toLong)
          } else {
            numOfWorkflowsToDeleteMetadataMetricActor ! CalculateNumberOfWorkflowsToDeleteMetadataMetricValue(currentTimestampMinusDelay)
          }
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

  def props(metadataDeletionConfig: MetadataDeletionConfig, serviceRegistryActor: ActorRef) = Props(new DeleteMetadataActor(metadataDeletionConfig, serviceRegistryActor))

  case object DeleteMetadataAction
}
