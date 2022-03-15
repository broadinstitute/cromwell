package cromwell.services.metadata.impl

import akka.actor.SupervisorStrategy.{Decider, Directive, Escalate, Resume}
import akka.actor.{Actor, ActorContext, ActorInitializationException, ActorLogging, ActorRef, Cancellable, OneForOneStrategy, Props}
import akka.routing.Listen
import cats.data.NonEmptyList
import com.typesafe.config.Config
import common.exception.AggregatedMessageException
import common.validation.Validation._
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.{LoadConfig, WorkflowId}
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.impl.MetadataDatabaseAccess.WorkflowArchiveStatusAndEndTimestamp
import cromwell.services.metadata.impl.MetadataStatisticsRecorder.MetadataStatisticsRecorderSettings
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import cromwell.services.metadata.impl.archiver.{ArchiveMetadataConfig, ArchiveMetadataSchedulerActor}
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.deleter.{DeleteMetadataActor, DeleteMetadataConfig}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object MetadataServiceActor {
  val MetadataInstrumentationPrefix = NonEmptyList.of("metadata")

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(MetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}

case class MetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends Actor
    with ActorLogging
    with MetadataDatabaseAccess
    with MetadataServicesStore
    with GracefulShutdownHelper
    with CromwellInstrumentation {

  private val decider: Decider = {
    case _: ActorInitializationException => Escalate
    case _ => Resume
  }

  override val supervisorStrategy = new OneForOneStrategy()(decider) {
    override def logFailure(context: ActorContext, child: ActorRef, cause: Throwable, decision: Directive) = {
      val childName = if (child == readActor) "Read" else "Write"
      log.error(cause, s"The $childName Metadata Actor died unexpectedly, metadata events might have been lost. Restarting it...")
    }
  }

  private val metadataSummaryRefreshInterval: Option[FiniteDuration] = {
    val duration = serviceConfig.getOrElse[Duration]("metadata-summary-refresh-interval", default = 1 second)
    if (duration.isFinite()) Option(duration.asInstanceOf[FiniteDuration]) else None
  }

  private val metadataSummaryRefreshLimit = serviceConfig.getOrElse("metadata-summary-refresh-limit", default = 5000)

  private val metadataReadTimeout: Duration =
    serviceConfig.getOrElse[Duration]("metadata-read-query-timeout", Duration.Inf)
  private val metadataReadRowNumberSafetyThreshold: Int =
    serviceConfig.getOrElse[Int]("metadata-read-row-number-safety-threshold", 1000000)

  private val metadataTableMetricsInterval: Option[FiniteDuration] = serviceConfig.getAs[FiniteDuration]("metadata-table-metrics-interval")

  private val metadataTableMetricsPath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "table"
  private val dataFreeMetricsPath: NonEmptyList[String] = metadataTableMetricsPath :+ "data_free"
  private val dataLengthMetricsPath: NonEmptyList[String] = metadataTableMetricsPath :+ "data_length"
  private val indexLengthMetricsPath:  NonEmptyList[String] = metadataTableMetricsPath :+ "index_length"

  def readMetadataWorkerActorProps(): Props =
    ReadDatabaseMetadataWorkerActor
      .props(metadataReadTimeout, metadataReadRowNumberSafetyThreshold)
      .withDispatcher(ServiceDispatcher)

  def metadataBuilderActorProps(): Props = MetadataBuilderActor
    .props(readMetadataWorkerActorProps, metadataReadRowNumberSafetyThreshold)
    .withDispatcher(ServiceDispatcher)

  val readActor = context.actorOf(ReadMetadataRegulatorActor.props(metadataBuilderActorProps, readMetadataWorkerActorProps), "ClassicMSA-ReadMetadataRegulatorActor")

  val dbFlushRate = serviceConfig.getOrElse("db-flush-rate", 5.seconds)
  val dbBatchSize = serviceConfig.getOrElse("db-batch-size", 200)
  val metadataWriteStatisticsConfig = MetadataStatisticsRecorderSettings(serviceConfig.as[Option[Config]]("metadata-write-statistics"))
  val writeActor = context.actorOf(WriteMetadataActor.props(dbBatchSize, dbFlushRate, serviceRegistryActor, LoadConfig.MetadataWriteThreshold, metadataWriteStatisticsConfig), "WriteMetadataActor")

  implicit val ec = context.dispatcher
  //noinspection ActorMutableStateInspection
  private var summaryRefreshCancellable: Option[Cancellable] = None

  private val summaryActor: Option[ActorRef] = buildSummaryActor

  summaryActor foreach { _ => self ! RefreshSummary }

  private val archiveMetadataActor: Option[ActorRef] = buildArchiveMetadataActor

  private val deleteMetadataActor: Option[ActorRef] = buildDeleteMetadataActor

  // if `metadata-table-size-metrics-interval` is specified, schedule sending size metrics at that interval
  metadataTableMetricsInterval.map(context.system.scheduler.schedule(1.minute, _, self, SendMetadataTableSizeMetrics)(context.dispatcher, self))

  private def scheduleSummary(): Unit = {
    metadataSummaryRefreshInterval foreach { interval =>
      summaryRefreshCancellable = Option(context.system.scheduler.scheduleOnce(interval, self, RefreshSummary)(context.dispatcher, self))
    }
  }

  override def postStop(): Unit = {
    summaryRefreshCancellable foreach { _.cancel() }
    super.postStop()
  }

  private def buildSummaryActor: Option[ActorRef] = {
    val actor = metadataSummaryRefreshInterval map {
      _ => context.actorOf(MetadataSummaryRefreshActor.props(serviceRegistryActor), "metadata-summary-actor")
    }
    val message = metadataSummaryRefreshInterval match {
      case Some(interval) => s"Metadata summary refreshing every $interval."
      case None => "Metadata summary refresh is off."
    }
    log.info(message)
    actor
  }

  private def buildArchiveMetadataActor: Option[ActorRef] = {
    if (serviceConfig.hasPath("archive-metadata")) {
      log.info("Building metadata archiver from config")
      ArchiveMetadataConfig.parseConfig(serviceConfig.getConfig("archive-metadata"))(context.system) match {
        case Right(config) => Option(context.actorOf(ArchiveMetadataSchedulerActor.props(config, serviceRegistryActor), "archive-metadata-scheduler"))
        case Left(errorList) => throw AggregatedMessageException("Failed to parse the archive-metadata config", errorList.toList)
      }
    } else {
      log.info("No metadata archiver defined in config")
      None
    }
  }

  private def buildDeleteMetadataActor: Option[ActorRef] = {
    if (serviceConfig.hasPath("delete-metadata")) {
      log.info("Building metadata deleter from config")
      DeleteMetadataConfig.parseConfig(serviceConfig.getConfig("delete-metadata")) match {
        case Right(config) => Option(context.actorOf(DeleteMetadataActor.props(config, serviceRegistryActor), "delete-metadata-actor"))
        case Left(errorList) => throw AggregatedMessageException("Failed to parse the archive-metadata config", errorList.toList)
      }
    } else {
      log.info("No metadata deleter defined in config")
      None
    }
  }

  private def validateWorkflowIdInMetadata(possibleWorkflowId: WorkflowId, sender: ActorRef): Unit = {
    workflowWithIdExistsInMetadata(possibleWorkflowId.toString) onComplete {
      case Success(true) => sender ! RecognizedWorkflowId
      case Success(false) => sender ! UnrecognizedWorkflowId
      case Failure(e) => sender ! FailedToCheckWorkflowId(new RuntimeException(s"Failed lookup attempt for workflow ID $possibleWorkflowId", e))
    }
  }

  private def validateWorkflowIdInMetadataSummaries(possibleWorkflowId: WorkflowId, sender: ActorRef): Unit = {
    workflowWithIdExistsInMetadataSummaries(possibleWorkflowId.toString) onComplete {
      case Success(true) => sender ! RecognizedWorkflowId
      case Success(false) => sender ! UnrecognizedWorkflowId
      case Failure(e) => sender ! FailedToCheckWorkflowId(new RuntimeException(s"Failed lookup attempt for workflow ID $possibleWorkflowId", e))
    }
  }

  private def fetchWorkflowMetadataArchiveStatusAndEndTime(workflowId: WorkflowId, sender: ActorRef): Unit = {
    getMetadataArchiveStatusAndEndTime(workflowId) onComplete {
      case Success(WorkflowArchiveStatusAndEndTimestamp(status, endTime)) =>
        MetadataArchiveStatus.fromDatabaseValue(status).toTry match {
          case Success(archiveStatus) => sender ! WorkflowMetadataArchivedStatusAndEndTime(archiveStatus, endTime)
          case Failure(e) => sender ! FailedToGetArchiveStatusAndEndTime(new RuntimeException(s"Failed to get metadata archive status for workflow ID $workflowId", e))
        }
      case Failure(e) => sender ! FailedToGetArchiveStatusAndEndTime(new RuntimeException(s"Failed to get metadata archive status for workflow ID $workflowId", e))
    }
  }

  private def sendMetadataTableSizeMetrics(): Unit = {
    getMetadataTableSizeInformation onComplete {
      case Success(v) =>
        v foreach { d =>
          sendGauge(dataLengthMetricsPath, d.dataLength)
          sendGauge(indexLengthMetricsPath, d.indexLength)
          sendGauge(dataFreeMetricsPath, d.dataFree)
        }
      case Failure(e) => log.error(e, s"Error fetching metadata table size metrics. Will try again in $metadataTableMetricsInterval...")
    }
  }

  def summarizerReceive: Receive = {
    case RefreshSummary => summaryActor foreach { _ ! SummarizeMetadata(metadataSummaryRefreshLimit, sender()) }
    case MetadataSummarySuccess => scheduleSummary()
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary()
  }

  def receive = summarizerReceive orElse {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(writeActor) ++ archiveMetadataActor.toList ++ deleteMetadataActor.toList)
    case SendMetadataTableSizeMetrics => sendMetadataTableSizeMetrics()
    case action: PutMetadataAction => writeActor forward action
    case action: PutMetadataActionAndRespond => writeActor forward action
    // Assume that listen messages are directed to the write metadata actor
    case listen: Listen => writeActor forward listen
    case v: ValidateWorkflowIdInMetadata => validateWorkflowIdInMetadata(v.possibleWorkflowId, sender())
    case v: ValidateWorkflowIdInMetadataSummaries => validateWorkflowIdInMetadataSummaries(v.possibleWorkflowId, sender())
    case g: FetchWorkflowMetadataArchiveStatusAndEndTime => fetchWorkflowMetadataArchiveStatusAndEndTime(g.workflowId, sender())
    case action: BuildMetadataJsonAction => readActor forward action
    case streamAction: GetMetadataStreamAction => readActor forward streamAction
  }
}
