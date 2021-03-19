package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef, FSM, LoggingFSM, Props}
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.Archived
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor.ArchiveWorkflowComplete
import cromwell.services.metadata.impl.archiver.StreamMetadataToGcsActor._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}


class StreamMetadataToGcsActor(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistry: ActorRef)
  extends Actor
    with ActorLogging
    with GracefulShutdownHelper
    with LoggingFSM[ArchivingMetadataState, ArchivingMetadataData]
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  implicit val ec: ExecutionContext = context.dispatcher

  startWith(Pending, NoData)

  when(Pending) {
    case Event(ArchiveMetadataForWorkflow(_), NoData) =>
      // TODO: Send a message to SRA to get a metadata stream
      ???
  }

  when(FetchingMetadata) {
    // TODO: Handle metadata streams by starting a asynchronous stream to GCS & pipe the results to ourselves in the
    // archiving workflow state below
    case _ => ???
  }

  when(ArchivingWorkflow) {
    case Event(MetadataArchiveResult(Success(_)), ArchivingData(workflowId)) =>
      scheduleDatabaseUpdateAndAwaitResult(workflowId, Archived)
    case Event(MetadataArchiveResult(Failure(_)), ArchivingData(_)) =>
      // TODO: Handle streaming failure
      ???
  }

  when(UpdatingDatabase) {
    case Event(DatabaseUpdateCompleted(Success(_)), UpdatingDatabaseData(workflowId, result)) =>
      // All set; send 'complete' to ArchiveMetadataSchedulerActor; wait for next workflow to archive
      context.parent ! ArchiveWorkflowComplete(workflowId, result)
      // Prepare for the next work request by reverting to the initial FSM state (and data):
      goto(Pending) using NoData
    case Event(DatabaseUpdateCompleted(Failure(reason)), UpdatingDatabaseData(workflowId, archivingResult)) =>
      log.error(reason, s"Unable to update workflow ID $workflowId's METADATA_ARCHIVE_STATUS to $archivingResult. Will try again in 1 minute")
      scheduleDatabaseUpdateAndAwaitResult(workflowId, archivingResult, delay = Option(1.minute))
  }

  whenUnhandled {
    case Event(ShutdownCommand, _) =>
      context stop self // TODO: Anything more graceful?
      stay()
    case other =>
      log.error(s"Programmer Error: Unexpected message to ${getClass.getSimpleName} ${self.path.name} in state $stateName with $stateData: $other")
      stay()
  }

  def scheduleDatabaseUpdateAndAwaitResult(workflowId: WorkflowId,
                                           newStatus: MetadataArchiveStatus,
                                           delay: Option[FiniteDuration] = None): FSM.State[ArchivingMetadataState, ArchivingMetadataData] = {

    def updateDatabase(): Unit = {
      val dbUpdateFuture = updateMetadataArchiveStatus(workflowId, newStatus)
      dbUpdateFuture onComplete { dbUpdateResult => self ! DatabaseUpdateCompleted(dbUpdateResult) }
    }

    delay match {
      case Some(d) => context.system.scheduler.scheduleOnce(d)(updateDatabase())
      case None => updateDatabase()
    }
    goto(UpdatingDatabase) using UpdatingDatabaseData(workflowId, newStatus)
  }
}

object StreamMetadataToGcsActor {
  def props(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistry: ActorRef): Props =
    Props(new StreamMetadataToGcsActor(archiveMetadataConfig, serviceRegistry))

  sealed trait StreamMetadataToGcsActorMessage
  final case class ArchiveMetadataForWorkflow(workflow: WorkflowId) extends StreamMetadataToGcsActorMessage
  final case class MetadataArchiveResult(archiveResult: Try[Unit]) extends StreamMetadataToGcsActorMessage
  final case class DatabaseUpdateCompleted(databaseUpdateResult: Try[Int]) extends StreamMetadataToGcsActorMessage

  sealed trait ArchivingMetadataState
  case object Pending extends ArchivingMetadataState
  case object FetchingMetadata extends ArchivingMetadataState
  case object ArchivingWorkflow extends ArchivingMetadataState
  case object UpdatingDatabase extends ArchivingMetadataState

  sealed trait ArchivingMetadataData
  case object NoData extends ArchivingMetadataData
  final case class FetchingData(workflowId: WorkflowId) extends ArchivingMetadataData
  final case class ArchivingData(workflowId: WorkflowId) extends ArchivingMetadataData
  final case class UpdatingDatabaseData(workflowId: WorkflowId, archiveResult: MetadataArchiveStatus) extends ArchivingMetadataData
}
