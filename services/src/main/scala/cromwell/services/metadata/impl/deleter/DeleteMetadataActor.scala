package cromwell.services.metadata.impl.deleter

import java.time.{OffsetDateTime, Duration => JDuration}
import java.util.concurrent.TimeUnit

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.util.Timeout
import cats.data.NonEmptyList
import common.util.StringUtil.EnhancedToStringable
import cromwell.core.instrumentation.InstrumentationPrefixes.ServicesPrefix
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.{Archived, ArchivedAndDeleted}
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}
import cromwell.services.metadata.impl.deleter.DeleteMetadataActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

class DeleteMetadataActor(deleteMetadataConfig: DeleteMetadataConfig,
                          override val serviceRegistryActor: ActorRef)
  extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher
  implicit val askTimeout: Timeout = new Timeout(60.seconds)

  private val deleterMetricsBasePath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "deleter"
  private val rowsDeletedMetricPath: NonEmptyList[String] = deleterMetricsBasePath :+ "rows_deleted"
  private val workflowsDeletedSuccessMetricPath: NonEmptyList[String] = deleterMetricsBasePath :+ "workflows_deleted" :+ "success"
  private val workflowsDeletedFailureMetricPath: NonEmptyList[String] = deleterMetricsBasePath :+ "workflows_deleted" :+ "failure"
  private val workflowDeleteTotalTimeMetricPath: NonEmptyList[String] = deleterMetricsBasePath :+ "workflow_delete_total_time"
  private val workflowsToDeleteMetricPath: NonEmptyList[String] = deleterMetricsBasePath :+ "workflows_to_delete"

  // Send an initial delete message to get things started:
  self ! DeleteNextWorkflowMessage

  // initial schedule for workflows left to delete metric
  context.system.scheduler.scheduleOnce(deleteMetadataConfig.instrumentationInterval)(workflowsLeftToDeleteMetric())

  override def receive: Receive = {
    case DeleteNextWorkflowMessage => {
      val startTime = OffsetDateTime.now()
      def calculateTimeSinceStart() = {
        FiniteDuration(JDuration.between(startTime, OffsetDateTime.now()).toMillis, TimeUnit.MILLISECONDS)
      }

      // These handlers send metrics for most paths even when they're not incremented, so that the metrics
      // paths are actively receiving data points throughout:
      deleteNextWorkflow().onComplete({
        case Success(true) =>
          increment(workflowsDeletedSuccessMetricPath, ServicesPrefix)
          count(workflowsDeletedFailureMetricPath, 0L, ServicesPrefix)
          sendTiming(workflowDeleteTotalTimeMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          self ! DeleteNextWorkflowMessage
        case Success(false) =>
          count(rowsDeletedMetricPath, 0L, ServicesPrefix)
          count(workflowsDeletedSuccessMetricPath, 0L, ServicesPrefix)
          count(workflowsDeletedFailureMetricPath, 0L, ServicesPrefix)
          sendGauge(workflowsToDeleteMetricPath, 0L, ServicesPrefix)
          sendTiming(workflowDeleteTotalTimeMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          scheduleNextDeleteAttemptAfterInterval()
          if (deleteMetadataConfig.debugLogging) log.info(s"No archived workflows which finished over ${deleteMetadataConfig.delayAfterWorkflowCompletion} ago remain to be deleted. Scheduling next poll in ${deleteMetadataConfig.backoffInterval}.")
        case Failure(error) =>
          count(rowsDeletedMetricPath, 0L, ServicesPrefix)
          count(workflowsDeletedSuccessMetricPath, 0L, ServicesPrefix)
          increment(workflowsDeletedFailureMetricPath, ServicesPrefix)
          sendTiming(workflowDeleteTotalTimeMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          log.error(error, s"Error while deleting, will wait ${deleteMetadataConfig.backoffInterval} then try again.")
          scheduleNextDeleteAttemptAfterInterval()
      })

    }
    case ShutdownCommand => context.stop(self)  // TODO: cancel any deletion action that might be happening?
    case other => log.info(s"Programmer Error! The DeleteMetadataSchedulerActor received unexpected message! ($sender sent ${other.toPrettyElidedString(1000)}})")
  }

  def workflowsLeftToDeleteMetric(): Unit = {
    val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(deleteMetadataConfig.delayAfterWorkflowCompletion.toSeconds)
    countWorkflowsLeftToDeleteThatEndedOnOrBeforeThresholdTimestamp(currentTimestampMinusDelay).onComplete({
      case Success(workflowsToDelete) =>
        sendGauge(workflowsToDeleteMetricPath, workflowsToDelete.longValue(), ServicesPrefix)
        // schedule next workflows left to delete query after interval
        context.system.scheduler.scheduleOnce(deleteMetadataConfig.instrumentationInterval)(workflowsLeftToDeleteMetric())
      case Failure(exception) =>
        log.error(exception, s"Something went wrong while fetching number of workflows left to delete. " +
          s"Scheduling next poll in ${deleteMetadataConfig.instrumentationInterval}.")
        // schedule next workflows left to delete query after interval
        context.system.scheduler.scheduleOnce(deleteMetadataConfig.instrumentationInterval)(workflowsLeftToDeleteMetric())
    })
  }

  def deleteNextWorkflow(): Future[Boolean] = for {
    maybeWorkflowId <- lookupNextWorkflowToDelete()
    result <- maybeWorkflowId match {
      case Some(id) =>
        log.info(s"Workflow $id identified for metadata deletion")
        for {
          rowsDeleted <- deleteAllMetadataEntriesForWorkflowAndUpdateArchiveStatus(id, MetadataArchiveStatus.toDatabaseValue(ArchivedAndDeleted))
          _ = count(rowsDeletedMetricPath, rowsDeleted.longValue(), ServicesPrefix)
          _ = log.info(s"Deleted $rowsDeleted metadata rows for $id")
        } yield true
      case None => Future.successful(false)
    }
  } yield result

  def lookupNextWorkflowToDelete(): Future[Option[WorkflowId]] = {
    val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(deleteMetadataConfig.delayAfterWorkflowCompletion.toSeconds)
    queryWorkflowIdsByArchiveStatusAndOlderThanTimestamp(
      MetadataArchiveStatus.toDatabaseValue(Archived),
      currentTimestampMinusDelay,
      batchSize = 1
    ).map(_.headOption.map(WorkflowId.fromString))
  }

  def scheduleNextDeleteAttemptAfterInterval(): Unit = {
    context.system.scheduler.scheduleOnce(deleteMetadataConfig.backoffInterval)(self ! DeleteNextWorkflowMessage)
    ()
  }

}

object DeleteMetadataActor {
  case object DeleteNextWorkflowMessage

  def props(config: DeleteMetadataConfig, serviceRegistryActor: ActorRef): Props = Props(new DeleteMetadataActor(config, serviceRegistryActor))
}
