package cromwell.services.metadata.hybridcarbonite

import java.nio.file.StandardOpenOption

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cromwell.core.WorkflowId
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataArchiveStatus.{ArchiveFailed, Archived}
import cromwell.services.metadata.MetadataService.GetMetadataAction
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowComplete
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor._
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.{SuccessfulMetadataJsonResponse, FailedMetadataJsonResponse}
import cromwell.services.metadata.{MetadataArchiveStatus, MetadataQuery}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

class CarbonitingMetadataFreezerActor(carboniterConfig: HybridCarboniteConfig,
                                      carboniteWorkerActor: ActorRef,
                                      serviceRegistry: ActorRef,
                                      ioActor: ActorRef) extends
  LoggingFSM[CarbonitingMetadataFreezingState, CarbonitingMetadataFreezingData]
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  implicit val ec: ExecutionContext = context.dispatcher

  val asyncIo = new AsyncIo(ioActor, DefaultIoCommandBuilder)

  startWith(Pending, NoData)

  when(Pending) {
    case Event(FreezeMetadata(workflowId), NoData) =>
      val fetchMetadataRequestToServiceRegistry = GetMetadataAction(MetadataQuery(
        workflowId = workflowId,
        jobKey = None,
        key = None,
        includeKeysOption = None,
        excludeKeysOption = Option(NonEmptyList.of("labels")),
        expandSubWorkflows = true
      ))

      serviceRegistry ! fetchMetadataRequestToServiceRegistry
      goto(Fetching) using FetchingData(workflowId)
  }

  when(Fetching) {
    case Event(SuccessfulMetadataJsonResponse(_, responseJson), FetchingData(workflowId)) =>
      asyncIo.writeAsync(carboniterConfig.makePath(workflowId), responseJson.prettyPrint, Seq(StandardOpenOption.CREATE)) onComplete {
        result => self ! CarbonitingFreezeResult(result)
      }
      goto(Freezing) using FreezingData(workflowId)

    case Event(FailedMetadataJsonResponse(_, reason), FetchingData(workflowId)) =>
      log.error(reason, s"Failed to fetch workflow $workflowId's metadata to archive. Marking as $ArchiveFailed")
      scheduleDatabaseUpdateAndAwaitResult(workflowId, ArchiveFailed)
  }

  when(Freezing) {
    case Event(CarbonitingFreezeResult(freezeResult), FreezingData(workflowId)) =>
      freezeResult.failed.foreach { reason => log.error(reason, s"Failed to freeze workflow $workflowId's metadata to GCS archive. Marking as $ArchiveFailed") }

      val newStatus = statusFromFreezeResult(freezeResult)
      scheduleDatabaseUpdateAndAwaitResult(workflowId, newStatus)
  }

  when(UpdatingDatabase) {
    case Event(DatabaseUpdateCompleted(Success(_)), UpdatingDatabaseData(workflowId, result)) =>
      // All set; send 'complete' to CarboniteWorkerActor; wait for next workflow to carbonite
      carboniteWorkerActor ! CarboniteWorkflowComplete(workflowId, result)
      // Prepare for the next work request by reverting to the initial FSM state (and data):
      goto(Pending) using NoData
    case Event(DatabaseUpdateCompleted(Failure(reason)), UpdatingDatabaseData(workflowId, freezingResult)) =>
      log.error(reason, s"Unable to update workflow ID $workflowId's METADATA_ARCHIVE_STATUS to $freezingResult. Will try again in 1 minute")
      scheduleDatabaseUpdateAndAwaitResult(workflowId, freezingResult, delay = Option(1.minute))
  }

  whenUnhandled {
    case Event(ShutdownCommand, _) =>
      context stop self
      stay()
    case other =>
      log.error(s"Programmer Error: Unexpected message to ${getClass.getSimpleName} ${self.path.name} in state $stateName with $stateData: $other")
      stay()
  }

  def statusFromFreezeResult(freezeResult: Try[Unit]): MetadataArchiveStatus = freezeResult match {
    case Success(_) => Archived
    case Failure(_) => ArchiveFailed
  }

  def scheduleDatabaseUpdateAndAwaitResult(workflowId: WorkflowId, newStatus: MetadataArchiveStatus, delay: Option[FiniteDuration] = None) = {

    def updateDatabase() = {
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


object CarbonitingMetadataFreezerActor {

  def props(carboniterConfig: HybridCarboniteConfig, carboniteWorkerActor: ActorRef, serviceRegistry: ActorRef, ioActor: ActorRef) =
    Props(new CarbonitingMetadataFreezerActor(carboniterConfig, carboniteWorkerActor, serviceRegistry, ioActor))

  sealed trait CarbonitingMetadataFreezingState
  case object Pending extends CarbonitingMetadataFreezingState
  case object Fetching extends CarbonitingMetadataFreezingState
  case object Freezing extends CarbonitingMetadataFreezingState
  case object UpdatingDatabase extends CarbonitingMetadataFreezingState

  sealed trait CarbonitingMetadataFreezingData
  case object NoData extends CarbonitingMetadataFreezingData
  final case class FetchingData(workflowId: WorkflowId) extends CarbonitingMetadataFreezingData
  final case class FreezingData(workflowId: WorkflowId) extends CarbonitingMetadataFreezingData
  final case class UpdatingDatabaseData(workflowId: WorkflowId, freezingResult: MetadataArchiveStatus) extends CarbonitingMetadataFreezingData

  sealed trait CarbonitingMetadataFreezingActorMessage
  final case class FreezeMetadata(workflow: WorkflowId) extends CarbonitingMetadataFreezingActorMessage
  final case class CarbonitingFreezeResult(freezingResult: Try[Unit]) extends CarbonitingMetadataFreezingActorMessage
  final case class DatabaseUpdateCompleted(databaseUpdateResult: Try[Int]) extends CarbonitingMetadataFreezingActorMessage

}
