package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}

import scala.util.{Failure, Success}
import scala.concurrent.duration._

class DeleteMetadataActor(metadataDeletionConfig: MetadataDeletionConfig, override val serviceRegistryActor: ActorRef) extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

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
          populateWorkflowsToDeleteMetadataMetric(currentTimestampMinusDelay)
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

  private def populateWorkflowsToDeleteMetadataMetric(currentTimestampMinusDelay: OffsetDateTime) = {
    // we need the dedicated `count` query here because the result of the `queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp`
    // query is upper-bounded by `batchSize` parameter, which is undesired here
    countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
      MetadataArchiveStatus.toDatabaseValue(Archived),
      currentTimestampMinusDelay
    ) onComplete {
      case Success(numToDelete) =>
        sendGauge(workflowsToDeleteMetadataMetricPath, numToDelete.toLong, InstrumentationPrefixes.ServicesPrefix)
      case Failure(ex) =>
        log.warning(s"Cannot send workflowsToDeleteMetadata gauge metric value: unable to query number of workflows-to-delete-metadata from the database: ${ex.getMessage}")
    }
  }
}

object DeleteMetadataActor {

  def props(metadataDeletionConfig: MetadataDeletionConfig, serviceRegistryActor: ActorRef) = Props(new DeleteMetadataActor(metadataDeletionConfig, serviceRegistryActor))

  case object DeleteMetadataAction

  private val workflowsToDeleteMetadataMetricPath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "delete" :+ "workflowsToDeleteMetadata"
}
