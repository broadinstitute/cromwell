package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.AsynchronousThrottlingGaugeMetricActor.{CalculateMetricValue, MetricValue}
import cromwell.services.instrumentation.{AsynchronousThrottlingGaugeMetricActor, CromwellInstrumentation}
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

class DeleteMetadataActor(metadataDeletionConfig: ActiveMetadataDeletionConfig, override val serviceRegistryActor: ActorRef) extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  log.info(s"Archived metadata deletion is configured to begin polling after ${metadataDeletionConfig.initialDelay}, and then delete up to ${metadataDeletionConfig.batchSize} workflows worth of metadata every ${metadataDeletionConfig.interval} (for workflows which completed at least ${metadataDeletionConfig.delayAfterWorkflowCompletion} ago).")
  context.system.scheduler.scheduleWithFixedDelay(metadataDeletionConfig.initialDelay, metadataDeletionConfig.interval, self, DeleteMetadataAction)

  private val workflowsToDeleteMetadataMetricPath =
    MetadataServiceActor.MetadataInstrumentationPrefix :+ "delete" :+ "numOfWorkflowsToDeleteMetadata"
  private val numOfWorkflowsToDeleteMetadataMetricActor =
    context.actorOf(AsynchronousThrottlingGaugeMetricActor.props(workflowsToDeleteMetadataMetricPath, InstrumentationPrefixes.ServicesPrefix, serviceRegistryActor))

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
          if (workflowIds.length < metadataDeletionConfig.batchSize) {
            numOfWorkflowsToDeleteMetadataMetricActor ! MetricValue(workflowIds.length.toLong)
          } else {
            def metricCalculatingFunction: ExecutionContext => Future[Int] =
              ec =>
                countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
                  MetadataArchiveStatus.toDatabaseValue(Archived),
                  currentTimestampMinusDelay
                )(ec)
            numOfWorkflowsToDeleteMetadataMetricActor ! CalculateMetricValue(metricCalculatingFunction)
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

  def props(metadataDeletionConfig: ActiveMetadataDeletionConfig, serviceRegistryActor: ActorRef) = Props(new DeleteMetadataActor(metadataDeletionConfig, serviceRegistryActor))
  case object DeleteMetadataAction
}
