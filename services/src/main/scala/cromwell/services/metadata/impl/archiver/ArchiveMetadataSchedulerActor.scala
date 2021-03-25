package cromwell.services.metadata.impl.archiver

import java.io.OutputStreamWriter
import java.nio.file.{Files, StandardOpenOption}

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.pattern.ask
import akka.util.Timeout
import common.util.StringUtil.EnhancedToStringable
import common.util.TimeUtil.EnhancedOffsetDateTime
import cromwell.core.path.Path
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.database.sql.SqlConverters.{ClobOptionToRawString, TimestampToSystemOffsetDateTime}
import cromwell.database.sql.tables.MetadataEntry
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataArchiveStatus.{Archived, Unarchived}
import cromwell.services.metadata.MetadataService.{GetMetadataStreamAction, MetadataLookupStreamFailed, MetadataLookupStreamSuccess, QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import slick.basic.DatabasePublisher

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

import java.time.OffsetDateTime


class ArchiveMetadataSchedulerActor(archiveMetadataConfig: ArchiveMetadataConfig,
                                    serviceRegistryActor: ActorRef)
  extends Actor
    with ActorLogging
    with GracefulShutdownHelper
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  implicit val ec: ExecutionContext = context.dispatcher
  implicit val askTimeout: Timeout = new Timeout(60.seconds)

  // kick off archiving immediately
  self ! ArchiveNextWorkflowMessage

  override def receive: Receive = {
    case ArchiveNextWorkflowMessage => archiveNextWorkflow().onComplete({
      case Success(true) => self ! ArchiveNextWorkflowMessage
      case Success(false) => scheduleNextWorkflowToArchive()
      case Failure(error) =>
        log.error(error, s"Error while archiving, will retry.")
        scheduleNextWorkflowToArchive()
    })
    case ShutdownCommand => context.stop(self)  // TODO: cancel any streaming that might be happening
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent $other})")
  }

  def archiveNextWorkflow(): Future[Boolean] = {
    for {
      maybeWorkflowId <- lookupNextWorkflowToArchive()
      result <- maybeWorkflowId match {
        case Some(id) => for {
          path <- Future(archiveMetadataConfig.makePath(id))
          _ = log.info(s"Archiving metadata for $id to ${path.pathAsString}")
          dbStream <- fetchStreamFromDatabase(id)
          _ <- streamMetadataToGcs(path, dbStream)
          _ <- updateMetadataArchiveStatus(id, Archived)
          _ = log.info(s"Archiving succeeded for $id")
        } yield true
        case None => Future.successful(false)
      }
    } yield result
  }

  def lookupNextWorkflowToArchive(): Future[Option[WorkflowId]] = {
    (serviceRegistryActor ? QueryForWorkflowsMatchingParameters(queryParametersForWorkflowsToArchive(OffsetDateTime.now(), archiveMetadataConfig.archiveDelay))) flatMap {
      case WorkflowQuerySuccess(response, _) =>
        if (response.results.nonEmpty)
          Future.successful(Option(WorkflowId.fromString(response.results.head.id)))
        else
          Future.successful(None)
      case WorkflowQueryFailure(reason) => Future.failed(new Exception("Failed to fetch new workflow to archive", reason))
      case other => Future.failed(new Exception(s"Programmer Error: Got unexpected message fetching new workflows to archive: ${other.toPrettyElidedString(1000)}"))
    }
  }

  def fetchStreamFromDatabase(workflowId: WorkflowId): Future[DatabasePublisher[MetadataEntry]] = {
    (serviceRegistryActor ? GetMetadataStreamAction(workflowId, archiveMetadataConfig.databaseStreamFetchSize)) flatMap {
      case MetadataLookupStreamSuccess(_, responseStream) => Future.successful(responseStream)
      case MetadataLookupStreamFailed(_, reason) => Future.failed(new Exception(s"Failed to get metadata stream", reason))
      case other => Future.failed(new Exception(s"Failed to get metadata stream: ${other.toPrettyElidedString(1000)}"))
    }
  }

  def streamMetadataToGcs(path: Path, stream: DatabasePublisher[MetadataEntry]): Future[Unit] = {
    for {
      csvPrinter <- Future {
        val gcsStream = Files.newOutputStream(path.nioPath, StandardOpenOption.CREATE)
        new CSVPrinter(new OutputStreamWriter(gcsStream), CSVFormat.DEFAULT.withHeader(CsvFileHeaders : _*))
      }
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

  def scheduleNextWorkflowToArchive(): Unit = {
    context.system.scheduler.scheduleOnce(archiveMetadataConfig.interval)(self ! ArchiveNextWorkflowMessage)
    ()
  }
}

object ArchiveMetadataSchedulerActor {
  case object ArchiveNextWorkflowMessage

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

  def queryParametersForWorkflowsToArchive(currentTime: OffsetDateTime, archiveDelay: FiniteDuration): Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "true",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1",
    NewestFirst.name -> "false", // oldest first for archiving
    EndDate.name -> currentTime.minusNanos(archiveDelay.toNanos).toUtcMilliString
  )

  def props(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistryActor: ActorRef): Props =
    Props(new ArchiveMetadataSchedulerActor(archiveMetadataConfig, serviceRegistryActor))
}
