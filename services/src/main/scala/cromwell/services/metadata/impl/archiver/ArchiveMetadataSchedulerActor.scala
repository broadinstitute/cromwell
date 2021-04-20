package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.pattern.ask
import akka.util.Timeout
import cats.data.NonEmptyList
import com.google.common.io.BaseEncoding
import com.google.common.primitives.Longs
import common.util.StringUtil.EnhancedToStringable
import common.util.TimeUtil.EnhancedOffsetDateTime
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.core.path.{Path, PathFactory}
import cromwell.core.instrumentation.InstrumentationPrefixes.ServicesPrefix
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.database.sql.SqlConverters.{ClobOptionToRawString, TimestampToSystemOffsetDateTime}
import cromwell.database.sql.tables.MetadataEntry
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus.{Archived, Unarchived}
import cromwell.services.metadata.MetadataService.{GetMetadataStreamAction, MetadataLookupStreamFailed, MetadataLookupStreamSuccess, QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQueryResult, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}
import cromwell.services.{IoActorRequester, MetadataServicesStore}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.codec.digest.PureJavaCrc32C
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import slick.basic.DatabasePublisher

import java.io.{OutputStream, OutputStreamWriter}
import java.nio.file.{Files, StandardOpenOption}
import java.time.{OffsetDateTime, Duration => JDuration}
import java.util.UUID
import java.util.concurrent.TimeUnit
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}


class ArchiveMetadataSchedulerActor(archiveMetadataConfig: ArchiveMetadataConfig,
                                    override val serviceRegistryActor: ActorRef)
  extends Actor
    with ActorLogging
    with GracefulShutdownHelper
    with MetadataDatabaseAccess
    with MetadataServicesStore
    with IoActorRequester
    with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher
  implicit val askTimeout: Timeout = new Timeout(60.seconds)
  lazy val futureAsyncIo: Future[AsyncIo] = requestIoActor() map { ioActor => {
    log.info(s"IoActor reference received by ${self.path.name}")
    new AsyncIo(ioActor, DefaultIoCommandBuilder)
  } }

  private val archiverMetricsBasePath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "archiver"
  private val rowsProcessedMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "rows_processed"
  private val workflowsProcessedSuccessMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_processed" :+ "success"
  private val workflowsProcessedFailureMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_processed" :+ "failure"
  private val workflowArchiveTotalTimeSuccessMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflow_archive_total_time" :+ "success"
  private val workflowArchiveTotalTimeFailureMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflow_archive_total_time" :+ "failure"
  private val workflowsToArchiveMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_to_archive"

  // kick off archiving immediately
  self ! ArchiveNextWorkflowMessage

  override def receive: Receive = {
    case ArchiveNextWorkflowMessage =>
      val startTime = OffsetDateTime.now()
      def calculateTimeSinceStart() = {
        FiniteDuration(JDuration.between(startTime, OffsetDateTime.now()).toMillis, TimeUnit.MILLISECONDS)
      }

      // These handlers send metrics for most paths even when they're not incremented, so that the metrics
      // paths are actively receiving data points throughout:
      archiveNextWorkflow().onComplete({
        case Success(true) =>
          increment(workflowsProcessedSuccessMetricPath, ServicesPrefix)
          count(workflowsProcessedFailureMetricPath, 0L, ServicesPrefix)
          sendTiming(workflowArchiveTotalTimeSuccessMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          self ! ArchiveNextWorkflowMessage
        case Success(false) =>
          count(rowsProcessedMetricPath, 0L, ServicesPrefix)
          count(workflowsProcessedSuccessMetricPath, 0L, ServicesPrefix)
          sendGauge(workflowsToArchiveMetricPath, 0L, ServicesPrefix)
          sendTiming(workflowArchiveTotalTimeSuccessMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          scheduleNextWorkflowToArchive()
          if (archiveMetadataConfig.debugLogging) log.info(s"No complete workflows which finished over ${archiveMetadataConfig.archiveDelay} ago remain to be archived. Scheduling next poll in ${archiveMetadataConfig.backoffInterval}.")
        case Failure(error) =>
          count(workflowsProcessedSuccessMetricPath, 0L, ServicesPrefix)
          increment(workflowsProcessedFailureMetricPath, ServicesPrefix)
          sendTiming(workflowArchiveTotalTimeFailureMetricPath, calculateTimeSinceStart(), ServicesPrefix)
          log.error(error, s"Error while archiving, will retry.")
          scheduleNextWorkflowToArchive()
      })
    case ShutdownCommand => context.stop(self)  // TODO: cancel any streaming that might be happening?
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent ${other.toPrettyElidedString(1000)}})")
  }

  def archiveNextWorkflow(): Future[Boolean] = {
    for {
      maybeWorkflowQueryResult <- lookupNextWorkflowToArchive()
      result <- maybeWorkflowQueryResult match {
        case Some(workflow) => for {
          path <- getGcsPathForMetadata(workflow)
          dbStream <- fetchStreamFromDatabase(WorkflowId(UUID.fromString(workflow.id)))
          _ = log.info(s"Archiving metadata for ${workflow.id} to ${path.pathAsString}")
          _ <- streamMetadataToGcs(path, dbStream)
          _ <- updateMetadataArchiveStatus(WorkflowId(UUID.fromString(workflow.id)), Archived)
          _ = log.info(s"Archiving succeeded for ${workflow.id}")
        } yield true
        case None => Future.successful(false)
      }
    } yield result
  }

  def lookupNextWorkflowToArchive(): Future[Option[WorkflowQueryResult]] = {
    (serviceRegistryActor ? QueryForWorkflowsMatchingParameters(queryParametersForWorkflowsToArchive(OffsetDateTime.now(), archiveMetadataConfig.archiveDelay))) flatMap {
      case WorkflowQuerySuccess(response, maybeQueryInfo) =>
        // Record total number of records available to archive (if provided in the response):
        for {
          queryInfo <- maybeQueryInfo
          totalRecords <- queryInfo.totalRecords
          _ = sendGauge(workflowsToArchiveMetricPath, totalRecords.longValue(), ServicesPrefix)
        } yield ()

        if (response.results.nonEmpty)
          Future.successful(Option(response.results.head))
        else
          Future.successful(None)
      case WorkflowQueryFailure(reason) => Future.failed(new Exception("Failed to fetch new workflow to archive", reason))
      case other => Future.failed(new Exception(s"Programmer Error: Got unexpected message fetching new workflows to archive: ${other.toPrettyElidedString(1000)}"))
    }
  }

  private def getGcsPathForMetadata(workflow: WorkflowQueryResult): Future[Path] =  {
    val bucket = archiveMetadataConfig.bucket
    val workflowId = workflow.id
    val rootWorkflowId = workflow.rootWorkflowId.getOrElse(workflowId)
    Future(PathFactory.buildPath(s"gs://$bucket/$rootWorkflowId/$workflowId.csv", archiveMetadataConfig.pathBuilders))
  }

  def fetchStreamFromDatabase(workflowId: WorkflowId): Future[DatabasePublisher[MetadataEntry]] = {
    (serviceRegistryActor ? GetMetadataStreamAction(workflowId)) flatMap {
      case MetadataLookupStreamSuccess(_, responseStream) => Future.successful(responseStream)
      case MetadataLookupStreamFailed(_, reason) => Future.failed(new Exception(s"Failed to get metadata stream", reason))
      case other => Future.failed(new Exception(s"Failed to get metadata stream: ${other.toPrettyElidedString(1000)}"))
    }
  }

  def streamMetadataToGcs(path: Path, stream: DatabasePublisher[MetadataEntry]): Future[Unit] = {
    val rowsCounter = new RowsCounterAndProgressiveLogger( logFunction = (newRows, totalRows) => {
      if (archiveMetadataConfig.debugLogging) logger.info(s"Uploaded $newRows new rows to ${path.pathAsString}. Total uploaded is now ${totalRows}") else ()
      count(rowsProcessedMetricPath, newRows, ServicesPrefix)
    }, 100000)
    for {
      asyncIo <- futureAsyncIo
      gcsStream = Files.newOutputStream(path.nioPath, StandardOpenOption.CREATE)
      crc32cStream = new Crc32cStream()
      teeStream = new TeeingOutputStream(gcsStream, crc32cStream)
      csvPrinter = new CSVPrinter(new OutputStreamWriter(teeStream), CSVFormat.DEFAULT.withHeader(CsvFileHeaders : _*))
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
        rowsCounter.increment()
      })
      _ = rowsCounter.finalLog()
      _ = csvPrinter.close()
      expectedChecksum = crc32cStream.checksumString
      uploadedChecksum <- asyncIo.hashAsync(path)
      _ <- if (uploadedChecksum == expectedChecksum) Future.successful(()) else Future.failed(new Exception(s"Uploaded checksum '$uploadedChecksum' did not match local calculation ('$expectedChecksum')"))
    } yield ()
  }

  def scheduleNextWorkflowToArchive(): Unit = {
    context.system.scheduler.scheduleOnce(archiveMetadataConfig.backoffInterval)(self ! ArchiveNextWorkflowMessage)
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

  final class TeeingOutputStream(streams: OutputStream*) extends OutputStream {
    override def write(b: Int): Unit = { streams.foreach(_.write(b)) }
    override def close(): Unit = { streams.foreach(_.close())}
    override def flush(): Unit = { streams.foreach(_.flush())}
  }

  final class Crc32cStream() extends OutputStream {
    private val checksumCalculator = new PureJavaCrc32C()
    override def write(b: Int): Unit = checksumCalculator.update(b)

    def checksumString: String = {
      val finalChecksumValue = checksumCalculator.getValue
      // Google checksums are actually only the lower four bytes of the crc32c checksum:
      val bArray = java.util.Arrays.copyOfRange(Longs.toByteArray(finalChecksumValue), 4, 8)
      BaseEncoding.base64.encode(bArray)
    }
  }

  final class RowsCounterAndProgressiveLogger(logFunction: (Long, Long) => Unit, logInterval: Int) {
    private var rowsSinceLog: Long = 0
    private var totalRows: Long = 0

    def increment(): Unit = {
      rowsSinceLog = rowsSinceLog + 1
      totalRows = totalRows + 1
      if (rowsSinceLog >= logInterval) {
        logFunction(rowsSinceLog, totalRows)
        rowsSinceLog = 0
      }
    }
    def finalLog(): Unit = {
      logFunction(rowsSinceLog, totalRows)
    }
  }
}
