package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus._
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor._
import cromwell.services.metadata.hybridcarbonite.WorkflowsToDeleteMetadataMetricHelperActor.{CalculateMetricAndSend, SendPrecalculatedMetric}
import cromwell.services.metadata.impl.MetadataDatabaseAccess

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

  val workflowsToDeleteMetadataMetricHelperActor =
    context.actorOf(WorkflowsToDeleteMetadataMetricHelperActor.props(self, serviceRegistryActor))

  override def receive: Receive = receive(isMetricActorFree = true)

  private def receive(isMetricActorFree: Boolean): Receive = {
    case MetricActorFreed => context.become(receive(isMetricActorFree = true))
    case DeleteMetadataAction =>
      val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(metadataDeletionConfig.delayAfterWorkflowCompletion.toSeconds)
      val workflowIdsForMetadataDeletionFuture = queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay,
        metadataDeletionConfig.batchSize
      )
      workflowIdsForMetadataDeletionFuture onComplete {
        case Success(workflowIds) =>
          if (isMetricActorFree) {
            log.info("WorkflowsToDeleteMetadataMetricHelperActor is free: sending a metric request.")
            context.become(receive(isMetricActorFree = false))
            if (workflowIds.length < metadataDeletionConfig.batchSize) {
              workflowsToDeleteMetadataMetricHelperActor ! SendPrecalculatedMetric(workflowIds.length.toLong)
            } else {
              workflowsToDeleteMetadataMetricHelperActor ! CalculateMetricAndSend(currentTimestampMinusDelay)
            }
          } else {
            log.info("WorkflowsToDeleteMetadataMetricHelperActor is busy: not sending another metric request at this time.")
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
  case object MetricActorFreed
}
