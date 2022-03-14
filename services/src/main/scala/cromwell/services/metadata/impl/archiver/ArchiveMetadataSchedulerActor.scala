package cromwell.services.metadata.impl.archiver

import java.io.{OutputStream, OutputStreamWriter}
import java.nio.file.{Files, StandardOpenOption}
import java.time.OffsetDateTime
import java.util.UUID

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
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus.Archived
import cromwell.services.metadata.MetadataService.{GetMetadataStreamAction, MetadataLookupStreamFailed, MetadataLookupStreamSuccess}
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}
import cromwell.services.{IoActorRequester, MetadataServicesStore}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.codec.digest.PureJavaCrc32C
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import slick.basic.DatabasePublisher

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}


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
  private val rowsPerWorkflowMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "rows_per_workflow"
  private val bytesProcessedMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "bytes_processed"
  private val bytesPerWorkflowMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "bytes_per_workflow"
  private val workflowsProcessedSuccessMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_processed" :+ "success"
  private val workflowsProcessedFailureMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_processed" :+ "failure"
  private val timeBehindExpectedDelayMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "time_behind_expected_delay"
  private val workflowArchiveTotalTimeMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflow_archive_total_time"
  private val workflowsToArchiveMetricPath: NonEmptyList[String] = archiverMetricsBasePath :+ "workflows_to_archive"
  private val archiverTimingMetricsBasePath: NonEmptyList[String] = archiverMetricsBasePath :+ "timings"
  private val archiverStreamTimingMetricsBasePath: NonEmptyList[String] = archiverTimingMetricsBasePath :+ "streaming"

  private val TerminalWorkflowStatuses: List[String] = List(WorkflowSucceeded, WorkflowAborted, WorkflowFailed).map(_.toString)

  // kick off archiving immediately
  self ! ArchiveNextWorkflowMessage

  // initial schedule for workflows left to archive metric
  context.system.scheduler.scheduleOnce(archiveMetadataConfig.instrumentationInterval)(workflowsLeftToArchiveMetric())

  override def receive: Receive = {
    case ArchiveNextWorkflowMessage =>
      val startTime = OffsetDateTime.now()

      // These handlers send metrics for most paths even when they're not incremented, so that the metrics
      // paths are actively receiving data points throughout:
      archiveNextWorkflows.onComplete({
        case Success(cnt) =>
          count(workflowsProcessedSuccessMetricPath, cnt, ServicesPrefix)
          count(workflowsProcessedFailureMetricPath, 0L, ServicesPrefix)

          if (cnt > 0) {
            if (archiveMetadataConfig.debugLogging) log.info(s"Successfully archived $cnt workflows.")
            self ! ArchiveNextWorkflowMessage
          } else {
            count(rowsProcessedMetricPath, 0L, ServicesPrefix)
            sendGauge(workflowsToArchiveMetricPath, 0L, ServicesPrefix)
            sendTiming(workflowArchiveTotalTimeMetricPath, calculateTimeSince(startTime), ServicesPrefix)
            if (archiveMetadataConfig.debugLogging) log.info(s"No complete workflows which finished over ${archiveMetadataConfig.archiveDelay} ago remain to be archived. Scheduling next poll in ${archiveMetadataConfig.backoffInterval}.")
            scheduleNextWorkflowToArchive()
          }
        case Failure(error) =>
          count(rowsProcessedMetricPath, 0L, ServicesPrefix)
          count(workflowsProcessedSuccessMetricPath, 0L, ServicesPrefix)
          count(workflowsProcessedFailureMetricPath, archiveMetadataConfig.batchSize, ServicesPrefix)
          sendTiming(workflowArchiveTotalTimeMetricPath, calculateTimeSince(startTime), ServicesPrefix)
          log.error(error, s"Error while archiving, will retry.")
          scheduleNextWorkflowToArchive()
      })
    case ShutdownCommand => context.stop(self)  // TODO: cancel any streaming that might be happening?
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent ${other.toPrettyElidedString(1000)}})")
  }

  def workflowsLeftToArchiveMetric(): Unit = {
    val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(archiveMetadataConfig.archiveDelay.toSeconds)
    countWorkflowsLeftToArchiveThatEndedOnOrBeforeThresholdTimestamp(
      TerminalWorkflowStatuses,
      currentTimestampMinusDelay
    ).onComplete({
      case Success(workflowsToArchive) =>
        sendGauge(workflowsToArchiveMetricPath, workflowsToArchive.longValue(), ServicesPrefix)
        // schedule next workflows left to archive query after interval
        context.system.scheduler.scheduleOnce(archiveMetadataConfig.instrumentationInterval)(workflowsLeftToArchiveMetric())
      case Failure(exception) =>
        log.error(exception, s"Something went wrong while fetching number of workflows left to archive. " +
          s"Scheduling next poll in ${archiveMetadataConfig.instrumentationInterval}.")
        // schedule next workflows left to archive query after interval
        context.system.scheduler.scheduleOnce(archiveMetadataConfig.instrumentationInterval)(workflowsLeftToArchiveMetric())
    })
  }

  def archiveNextWorkflows: Future[Long] = {
    val batchLookupStartTime = OffsetDateTime.now()

    for {
      workflowSummaryEntries <- lookupNextWorkflowsToArchive(archiveMetadataConfig.batchSize)
      batchLookupEndTime = OffsetDateTime.now()
      _ = sendTiming(archiverTimingMetricsBasePath :+ "lookup_next_workflows", calculateTimeDifference(batchLookupStartTime, batchLookupEndTime), ServicesPrefix)
      _ = if (archiveMetadataConfig.debugLogging) log.info(s"About to archive batch of ${workflowSummaryEntries.size} workflows.")
      result <- archiveSummaryEntries(workflowSummaryEntries)
      _ = sendTiming(archiverTimingMetricsBasePath :+ "batch_lookup_and_archive_time", calculateTimeDifference(batchLookupStartTime, OffsetDateTime.now()), ServicesPrefix)
    } yield result
  }

  private def archiveSummaryEntries(entries: Seq[WorkflowMetadataSummaryEntry]): Future[Long] = {
    if (entries.isEmpty) {
      sendGauge(timeBehindExpectedDelayMetricPath, 0L, ServicesPrefix)
      Future.successful(0L)
    } else {
      val resultSeq: Seq[Future[Long]] = entries.map(archiveSummaryEntry)
      val result: Future[Seq[Long]] = Future.sequence(resultSeq)

      result.map(_.sum)
    }
  }

  private def archiveSummaryEntry(entry: WorkflowMetadataSummaryEntry): Future[Long] = {
    entry.endTimestamp.foreach { workflowEndTime =>
      sendGauge(timeBehindExpectedDelayMetricPath, calculateTimeSince(workflowEndTime.toSystemOffsetDateTime).toMillis - archiveMetadataConfig.archiveDelay.toMillis, ServicesPrefix)
    }

    val archiveStartTime = OffsetDateTime.now()

    for {
      path <- Future.fromTry(getGcsPathForMetadata(entry))
      workflowId = entry.workflowExecutionUuid
      dbStream <- fetchStreamFromDatabase(WorkflowId(UUID.fromString(workflowId)))
      _ = log.info(s"Archiving metadata for $workflowId to ${path.pathAsString}")
      readyToStreamTime = OffsetDateTime.now()
      _ = sendTiming(archiverTimingMetricsBasePath :+ "prepare_to_stream", calculateTimeDifference(archiveStartTime, readyToStreamTime), ServicesPrefix)
      _ <- streamMetadataToGcs(path, dbStream)
      streamCompleteTime = OffsetDateTime.now()
      _ = sendTiming(archiverTimingMetricsBasePath :+ "stream_to_gcs", calculateTimeDifference(readyToStreamTime, streamCompleteTime), ServicesPrefix)
      _ <- updateMetadataArchiveStatus(WorkflowId(UUID.fromString(workflowId)), Archived)
      statusUpdatedTime = OffsetDateTime.now()
      _ = sendTiming(archiverTimingMetricsBasePath :+ "archive_status_update", calculateTimeDifference(streamCompleteTime, statusUpdatedTime), ServicesPrefix)
      _ = sendTiming(workflowArchiveTotalTimeMetricPath, calculateTimeDifference(archiveStartTime, statusUpdatedTime), ServicesPrefix)
      _ = log.info(s"Archiving succeeded for $workflowId")
    } yield 1L
  }

  def lookupNextWorkflowsToArchive(count: Long): Future[Seq[WorkflowMetadataSummaryEntry]] = {
    val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(archiveMetadataConfig.archiveDelay.toSeconds)
    queryWorkflowsToArchiveThatEndedOnOrBeforeThresholdTimestamp(
      TerminalWorkflowStatuses,
      currentTimestampMinusDelay,
      batchSize = count
    )
  }

  private def getGcsPathForMetadata(summaryEntry: WorkflowMetadataSummaryEntry): Try[Path] =  {
    /*
      Note: The naming convention for archived workflows is:
        - if a workflow has no root workflow, its archived metadata is put in GCS at <workflow_id>/<workflow_id>.csv
        - if a workflow is a subworkflow, its archived metadata is put in GCS under it's root workflow's directory i.e.
          <root_workflow_id>/<subworkflow_id>.csv
      Changing this convention would break the expectations of where to find the archived metadata files.
   */
    val bucket = archiveMetadataConfig.bucket
    val workflowId = summaryEntry.workflowExecutionUuid
    val rootWorkflowId = summaryEntry.rootWorkflowExecutionUuid.getOrElse(workflowId)
    Try {
      PathFactory.buildPath(s"gs://$bucket/$rootWorkflowId/$workflowId.csv", archiveMetadataConfig.pathBuilders)
    }
  }

  def fetchStreamFromDatabase(workflowId: WorkflowId): Future[DatabasePublisher[MetadataEntry]] = {
    (serviceRegistryActor ? GetMetadataStreamAction(workflowId)) flatMap {
      case MetadataLookupStreamSuccess(_, responseStream) => Future.successful(responseStream)
      case MetadataLookupStreamFailed(_, reason) => Future.failed(new Exception(s"Failed to get metadata stream", reason))
      case other => Future.failed(new Exception(s"Failed to get metadata stream: ${other.toPrettyElidedString(1000)}"))
    }
  }

  def streamMetadataToGcs(path: Path, stream: DatabasePublisher[MetadataEntry]): Future[Unit] = {
    val streamStartTime = OffsetDateTime.now()

    val rowsCounter = new CounterAndProgressiveLogger( logFunction = (newRows, totalRows) => {
      if (archiveMetadataConfig.debugLogging) logger.info(s"Uploaded $newRows new rows to ${path.pathAsString}. Total uploaded is now ${totalRows}") else ()
      count(rowsProcessedMetricPath, newRows, ServicesPrefix)
    }, 100000)

    for {
      asyncIo <- futureAsyncIo
      gotAsyncIoTime = OffsetDateTime.now()
      _ = sendTiming(archiverStreamTimingMetricsBasePath :+ "get_async_io", calculateTimeDifference(streamStartTime, gotAsyncIoTime), ServicesPrefix)
      gcsStream = Files.newOutputStream(path.nioPath, StandardOpenOption.CREATE)
      gcsStreamCreatedTime = OffsetDateTime.now()
      _ = sendTiming(archiverStreamTimingMetricsBasePath :+ "create_gcs_stream", calculateTimeDifference(gotAsyncIoTime, gcsStreamCreatedTime), ServicesPrefix)
      crc32cStream = new Crc32cStream()
      teeStream = new TeeingOutputStream(gcsStream, crc32cStream, new ByteCountingOutputStream())
      csvPrinter =
        new CSVPrinter(
          new OutputStreamWriter(teeStream),
          CSVFormat.DEFAULT.builder().setHeader(CsvFileHeaders : _*).build(),
        )
      csvPrinterCreatedTime = OffsetDateTime.now()
      _ = sendTiming(archiverStreamTimingMetricsBasePath :+ "create_csv_printer", calculateTimeDifference(gcsStreamCreatedTime, csvPrinterCreatedTime), ServicesPrefix)
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
      _ = rowsCounter.manualLog()
      _ = sendGauge(rowsPerWorkflowMetricPath, rowsCounter.getTotalCount, ServicesPrefix)
      _ = csvPrinter.close()
      streamingCompleteTime = OffsetDateTime.now()
      _ = sendTiming(archiverStreamTimingMetricsBasePath :+ "stream_data_to_gcs", calculateTimeDifference(csvPrinterCreatedTime, streamingCompleteTime), ServicesPrefix)
      expectedChecksum = crc32cStream.checksumString
      uploadedChecksum <- asyncIo.hashAsync(path)
      checksumValidatedTime = OffsetDateTime.now()
      _ = sendTiming(archiverStreamTimingMetricsBasePath :+ "checksum_validation", calculateTimeDifference(streamingCompleteTime, checksumValidatedTime), ServicesPrefix)
      _ <- if (uploadedChecksum == expectedChecksum) Future.successful(()) else Future.failed(new Exception(s"Uploaded checksum '$uploadedChecksum' did not match local calculation ('$expectedChecksum')"))
    } yield ()
  }

  def scheduleNextWorkflowToArchive(): Unit = {
    context.system.scheduler.scheduleOnce(archiveMetadataConfig.backoffInterval)(self ! ArchiveNextWorkflowMessage)
    ()
  }

  final class ByteCountingOutputStream() extends OutputStream {
    val byteCounter = new CounterAndProgressiveLogger( logFunction = (newBytes, totalBytes) => {
      if (archiveMetadataConfig.debugLogging) logger.info(s"Uploaded $newBytes new bytes. Total uploaded is now $totalBytes") else ()
      count(bytesProcessedMetricPath, newBytes, ServicesPrefix)
    }, 100000)

    override def write(b: Int): Unit = byteCounter.increment()
    override def close(): Unit = {
      byteCounter.manualLog()
      sendGauge(bytesPerWorkflowMetricPath, byteCounter.getTotalCount, ServicesPrefix)
    }
    override def flush(): Unit = byteCounter.manualLog()
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

  final class CounterAndProgressiveLogger(logFunction: (Long, Long) => Unit, logInterval: Int) {
    private var countSinceLog: Long = 0
    private var totalCount: Long = 0

    def increment(): Unit = {
      countSinceLog = countSinceLog + 1
      totalCount = totalCount + 1
      if (countSinceLog >= logInterval) {
        logFunction(countSinceLog, totalCount)
        countSinceLog = 0
      }
    }
    def manualLog(): Unit = {
      logFunction(countSinceLog, totalCount)
      countSinceLog = 0
    }
    def getTotalCount: Long = totalCount
  }
}
