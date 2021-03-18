package cromwell.services.metadata.hybridcarbonite

import java.io.{OutputStream, OutputStreamWriter}
import java.nio.file.{Files, StandardOpenOption}

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import common.util.TimeUtil.EnhancedOffsetDateTime
import cromwell.core.WorkflowId
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.core.path.Path
import cromwell.database.sql.SqlConverters.{ClobOptionToRawString, TimestampToSystemOffsetDateTime}
import cromwell.database.sql.tables.MetadataEntry
import cromwell.services.metadata.MetadataArchiveStatus.{ArchiveFailed, Archived, TooLargeToArchive}
import cromwell.services.metadata.MetadataService.{GetMetadataStreamAction, MetadataLookupStreamResponse}
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowComplete
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor._
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.{MetadataArchiveStatus, MetadataQuery}
import cromwell.services.{FailedMetadataJsonResponse, MetadataServicesStore, MetadataTooLargeException}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import slick.basic.DatabasePublisher

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

class CarbonitingMetadataFreezerActor(freezingConfig: ActiveMetadataFreezingConfig,
                                      carboniterConfig: HybridCarboniteConfig,
                                      carboniteWorkerActor: ActorRef,
                                      serviceRegistry: ActorRef,
                                      ioActor: ActorRef) extends
  LoggingFSM[CarbonitingMetadataFreezingState, CarbonitingMetadataFreezingData]
  with MetadataDatabaseAccess
  with MetadataServicesStore {

  implicit val ec: ExecutionContext = context.dispatcher

  val asyncIo = new AsyncIo(ioActor, DefaultIoCommandBuilder)

  val CsvFileHeaders = List(
    "METADATA_JOURNAL_ID",
    "WORKFLOW_EXECUTION_UUID",
    "METADATA_KEY",
    "CALL_FQN",
    "JOB_SCATTER_INDEX",
    "JOB_RETRY_ATTEMPT",
    "METADATA_VALUE",
    "METADATA_TIMESTAMP",
    "METADATA_VALUE_TYPE"
  )

  startWith(Pending, NoData)

  when(Pending) {
    case Event(FreezeMetadata(workflowId), NoData) =>
      log.info(s"Submitting request for metadata stream for $workflowId")
      val fetchMetadataRequestToServiceRegistry = GetMetadataStreamAction(MetadataQuery(
        workflowId = workflowId,
        jobKey = None,
        key = None,
        includeKeysOption = None,
        excludeKeysOption = None,
        expandSubWorkflows = true
      ))

      serviceRegistry ! fetchMetadataRequestToServiceRegistry
      goto(Fetching) using FetchingData(workflowId)
  }

//  final case class CsvFriendlyMetadataEntry(
//    workflowExecutionUuid: String,
//    callFullyQualifiedName: String,
//    jobIndex: Int,
//    jobAttempt: Int,
//    metadataKey: String,
//    metadataValue: String,
//    metadataValueType: String,
//    metadataTimestamp: String
//  )
//  object CsvFriendlyMetadataEntry {
//    def apply(me: MetadataEntry): CsvFriendlyMetadataEntry = {
//      // Timestamp format YYYY-MM-DDTHH:MM:SS.00Z as per https://stackoverflow.com/questions/47466296/bigquery-datetime-format-csv-to-bigquery-yyyy-mm-dd-hhmmss-ssssss
//
//      val declobberedValue = me.metadataValue.map(_.toString).getOrElse("")
//      val bqFriendlyTimestamp = me.metadataTimestamp.formatted("YYYY-MM-DDTHH:MM:SS.00Z")
//
//      new CsvFriendlyMetadataEntry(
//        me.workflowExecutionUuid,
//        me.callFullyQualifiedName.getOrElse(""),
//        me.jobIndex.getOrElse(-1),
//        me.jobAttempt.getOrElse(1),
//        me.metadataKey,
//        declobberedValue, //me.metadataValue
//        me.metadataValueType.getOrElse(""),
//        bqFriendlyTimestamp
//      )
//    }
//  }

  def writeStreamToGcs(workflowId: WorkflowId, stream: DatabasePublisher[MetadataEntry]): Future[Unit] = {
    val getCsvPrinter: Future[CSVPrinter] = {
      Future {
        val path: Path = carboniterConfig.makePath(workflowId)
        val gcsStream: OutputStream = Files.newOutputStream(path.nioPath, StandardOpenOption.CREATE)
        val printWriter = new OutputStreamWriter(gcsStream)
        new CSVPrinter(printWriter, CSVFormat.DEFAULT.withHeader(CsvFileHeaders : _*))
      }
    }

    for {
      csvPrinter <- getCsvPrinter
      _ <- stream.foreach(me => {
        csvPrinter.printRecord(
          me.metadataEntryId.map(_.toString).getOrElse(""),
          me.workflowExecutionUuid,
          me.metadataKey,
          me.callFullyQualifiedName.getOrElse(""),
          me.jobIndex.map(_.toString).getOrElse(""),
          me.jobAttempt.map(_.toString).getOrElse(""),
          me.metadataValue.toRawString,
          me.metadataTimestamp.toSystemOffsetDateTime.toUtcMilliString,
          me.metadataValueType.getOrElse("")
        )
      })
      _ = csvPrinter.close()
    } yield ()
  }

  when(Fetching) {
    case Event(MetadataLookupStreamResponse(_, responseStream), FetchingData(workflowId)) =>
      log.info(s"Received metadata stream for $workflowId. Beginning stream...")

      writeStreamToGcs(workflowId, responseStream) onComplete {
        result => self ! CarbonitingFreezeResult(result)
      }
      goto(Freezing) using FreezingData(workflowId)

    case Event(FailedMetadataJsonResponse(_, reason: MetadataTooLargeException), FetchingData(workflowId)) =>
      log.error(reason, s"Carboniting failure: $reason. Marking as $TooLargeToArchive")
      scheduleDatabaseUpdateAndAwaitResult(workflowId, TooLargeToArchive)

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

  def scheduleDatabaseUpdateAndAwaitResult(workflowId: WorkflowId,
                                           newStatus: MetadataArchiveStatus,
                                           delay: Option[FiniteDuration] = None): FSM.State[CarbonitingMetadataFreezingState, CarbonitingMetadataFreezingData] = {

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


object CarbonitingMetadataFreezerActor {

  def props(freezingConfig: ActiveMetadataFreezingConfig, carboniterConfig: HybridCarboniteConfig, carboniteWorkerActor: ActorRef, serviceRegistry: ActorRef, ioActor: ActorRef) =
    Props(new CarbonitingMetadataFreezerActor(freezingConfig, carboniterConfig, carboniteWorkerActor, serviceRegistry, ioActor))

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
