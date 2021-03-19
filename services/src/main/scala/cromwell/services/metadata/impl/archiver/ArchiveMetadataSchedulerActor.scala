package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus.Unarchived
import cromwell.services.metadata.MetadataService.{QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor._
import cromwell.services.metadata.impl.archiver.StreamMetadataToGcsActor.ArchiveMetadataForWorkflow
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.{Duration, MILLISECONDS}
import scala.util.{Failure, Success, Try}


class ArchiveMetadataSchedulerActor(archiveMetadataConfig: ArchiveMetadataConfig,
                                    override val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher
  val streamMetadataActor: ActorRef = context.actorOf(StreamMetadataToGcsActor.props(archiveMetadataConfig, serviceRegistryActor))

  // Schedule the initial workflows to archive query
  scheduleNextWorkflowsToArchiveQuery()

  override def receive: Receive = {
    case querySuccess: WorkflowQuerySuccess =>
      // TODO: Instrument "workflows to archive" query time?

      if (querySuccess.response.results.nonEmpty)
        archiveWorkflow(querySuccess.response.results.head.id)
      else
        scheduleNextWorkflowsToArchiveQuery()
    case queryFailure: WorkflowQueryFailure =>
      log.error(queryFailure.reason, s"Error while querying workflow to archive, will retry.")
      scheduleNextWorkflowsToArchiveQuery()
    case archiveComplete: ArchiveWorkflowComplete =>
      // TODO: Instrument time to archive workflow?
      if (archiveMetadataConfig.debugLogging) { log.info(s"Archiving complete for workflow ${archiveComplete.workflowId}") }

      // Immediately reset the timer and check for the next workflow to archive:
      nothingToArchiveBackoff.googleBackoff.reset()
      queryForWorkflowToArchive()
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(streamMetadataActor))
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent $other})")
  }

  def archiveWorkflow(workflowId: String): Unit = {
    Try(WorkflowId.fromString(workflowId)) match {
      case Success(id: WorkflowId) =>
        if (archiveMetadataConfig.debugLogging) { log.info(s"Starting to archive workflow: $workflowId") }
        streamMetadataActor ! ArchiveMetadataForWorkflow(id)
      case Failure(e) =>
        log.error(e, s"Programmer Error: Cannot archive workflow $workflowId. Error while converting ${workflowId} to WorkflowId, will retry.")
        scheduleNextWorkflowsToArchiveQuery()
    }
  }

  val nothingToArchiveBackoff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = archiveMetadataConfig.initialInterval,
    maxInterval = archiveMetadataConfig.maxInterval,
    multiplier = archiveMetadataConfig.multiplier,
    randomizationFactor = 0.0
  )

  def scheduleNextWorkflowsToArchiveQuery(): Unit = {
    val duration = Duration(nothingToArchiveBackoff.backoffMillis, MILLISECONDS)
    context.system.scheduler.scheduleOnce(duration)(queryForWorkflowToArchive())
    ()
  }

  def queryForWorkflowToArchive(): Unit = {
    serviceRegistryActor ! QueryForWorkflowsMatchingParameters(queryParametersForWorkflowsToArchive)
  }
}

object ArchiveMetadataSchedulerActor {

  sealed trait ArchiveWorkflowWorkerMessage
  case class ArchiveWorkflowComplete(workflowId: WorkflowId, result: cromwell.services.metadata.MetadataArchiveStatus) extends ArchiveWorkflowWorkerMessage

  // TODO: Archive from oldest-first
  // TODO: Allow requirements like "End timestamp not within 1y (eg)"
  val queryParametersForWorkflowsToArchive: Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "true",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  )

  def props(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistryActor: ActorRef): Props =
    Props(new ArchiveMetadataSchedulerActor(archiveMetadataConfig, serviceRegistryActor))
}
