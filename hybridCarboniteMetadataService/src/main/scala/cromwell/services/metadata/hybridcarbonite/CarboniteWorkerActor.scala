package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus.{ArchiveFailed, Archived, Unarchived}
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{DeleteMetadataAction, DeleteMetadataFailedResponse, DeleteMetadataSuccessfulResponse, QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor._
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor.FreezeMetadata
import cromwell.services.metadata.impl.MetadataServiceActor
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

class CarboniteWorkerActor(carboniterConfig: HybridCarboniteConfig,
                           override val serviceRegistryActor: ActorRef,
                           ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher

  val carboniteFreezerActor = context.actorOf(CarbonitingMetadataFreezerActor.props(carboniterConfig, self, serviceRegistryActor, ioActor))

  val backOff = {
    val scanConfig = carboniterConfig.freezeScanConfig
    SimpleExponentialBackoff(
      initialInterval = scanConfig.initialInterval,
      maxInterval = scanConfig.maxInterval,
      multiplier = scanConfig.multiplier,
      randomizationFactor = 0.0
    )
  }

  scheduleNextCarbonitingQuery()

  var workflowQueryStartTime: Option[Long] = None
  var carboniteFreezeStartTime: Option[Long] = None

  override def receive: Receive = {
    case s: WorkflowQuerySuccess =>

      s.meta match {
        case Some(MetadataService.QueryMetadata(_, _, Some(totalRecords))) => sendGauge(CarboniteWorkflowsToCarboniteMetricPath, totalRecords.toLong, InstrumentationPrefix)
        case _ => log.warning("Unable to write a 'workflowsToCarbonite' metric. No 'totalRecords' count was returned by the workflows query")
      }

      recordTimeSinceStartAsMetric(CarboniteWorkflowQueryTimeMetricPath, workflowQueryStartTime)
      workflowQueryStartTime = None

      if (s.response.results.nonEmpty)
        carboniteWorkflow(s.response.results.head.id)
      else
        scheduleNextCarbonitingQuery()

    case f: WorkflowQueryFailure =>
      log.error(f.reason, s"Error while querying workflow to carbonite, will retry.")
      scheduleNextCarbonitingQuery()
    case c: CarboniteWorkflowComplete =>

      recordTimeSinceStartAsMetric(CarboniteFreezingTimeMetricPath, carboniteFreezeStartTime)
      carboniteFreezeStartTime = None

      if (carboniterConfig.debugLogging) { log.info(s"Carboniting complete for workflow ${c.workflowId}") }

      c.result match {
        case Archived =>
          increment(CarboniteSuccessesMetricPath, InstrumentationPrefix)

          if (carboniterConfig.debugLogging) { log.info(s"Starting deleting metadata from database for carbonited workflow: ${c.workflowId}") }
          serviceRegistryActor ! DeleteMetadataAction(c.workflowId, self)
        case ArchiveFailed =>
          increment(CarboniteFailuresMetricPath, InstrumentationPrefix)
        case _ =>
          log.error(s"Programmer Error! The CarboniteWorkerActor cannot convert this into a completion metric: $c")
      }

      if (c.result != Archived) {
        resetBackoffAndFindNextWorkflowForCarboniting()
      }
    case DeleteMetadataSuccessfulResponse(workflowId) =>
      if (carboniterConfig.debugLogging) { log.info(s"Completed deleting metadata from database for carbonited workflow: $workflowId") }
      resetBackoffAndFindNextWorkflowForCarboniting()
    case DeleteMetadataFailedResponse(workflowId, reason) =>
      log.error(reason, s"All attempts to delete metadata from database for carbonited workflow $workflowId failed.")
      resetBackoffAndFindNextWorkflowForCarboniting()
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(carboniteFreezerActor))
    case other => log.error(s"Programmer Error! The CarboniteWorkerActor received unexpected message! ($sender sent $other})")
  }

  private def resetBackoffAndFindNextWorkflowForCarboniting() = {
    // Immediately reset the timer and check for the next workflow to carbonite:
    backOff.googleBackoff.reset()
    findWorkflowToCarbonite()
  }

  def scheduleNextCarbonitingQuery(resetBackoff: Boolean = false): Unit = {
    if (resetBackoff) backOff.googleBackoff.reset()
    val duration = Duration(backOff.backoffMillis, MILLISECONDS)
    context.system.scheduler.scheduleOnce(duration)(findWorkflowToCarbonite())
    ()
  }


  def findWorkflowToCarbonite(): Unit = {
    workflowQueryStartTime = Option(System.currentTimeMillis())
    serviceRegistryActor ! QueryForWorkflowsMatchingParameters(CarboniteWorkerActor.buildQueryParametersForWorkflowToCarboniteQuery(carboniterConfig.minimumSummaryEntryId))
  }


  def carboniteWorkflow(workflowId: String): Unit = {
    carboniteFreezeStartTime = Option(System.currentTimeMillis())

    Try(WorkflowId.fromString(workflowId)) match {
      case Success(id: WorkflowId) =>
        if (carboniterConfig.debugLogging) { log.info(s"Starting carboniting of workflow: $workflowId") }
        carboniteFreezerActor ! FreezeMetadata(id)
      case Failure(e) =>
        log.error(e, s"Cannot carbonite workflow $workflowId. Error while converting it to WorkflowId, will retry.")
        scheduleNextCarbonitingQuery()
    }
  }

  def recordTimeSinceStartAsMetric(path: NonEmptyList[String], startTime: Option[Long]) = startTime match {
    case Some(timestamp) =>
      val duration = (System.currentTimeMillis() - timestamp).millis
      sendTiming(path, duration, InstrumentationPrefix)
    case None =>
      log.error(s"Programmer Error: CarboniteWorkerActor cannot record '${path.toList.mkString(":")}' metric. Start time was not set.")
  }
}

object CarboniteWorkerActor {

  sealed trait CarboniteWorkerMessage
  case class CarboniteWorkflowComplete(workflowId: WorkflowId, result: cromwell.services.metadata.MetadataArchiveStatus) extends CarboniteWorkerMessage

  def buildQueryParametersForWorkflowToCarboniteQuery(minimumSummaryEntryId: Option[Long]): Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "false",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  ) ++ minimumSummaryEntryId.map(id => MinimumSummaryEntryId.name -> s"$id")

  def props(carboniterConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef, ioActor: ActorRef) =
    Props(new CarboniteWorkerActor(carboniterConfig, serviceRegistryActor, ioActor))

  val InstrumentationPrefix: Option[String] = InstrumentationPrefixes.ServicesPrefix
  private val carboniteMetricsBasePath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "carbonite"

  val CarboniteSuccessesMetricPath: NonEmptyList[String] = carboniteMetricsBasePath :+ "successful"
  val CarboniteFailuresMetricPath: NonEmptyList[String] = carboniteMetricsBasePath :+ "failure"

  val CarboniteWorkflowQueryTimeMetricPath: NonEmptyList[String] = carboniteMetricsBasePath :+ "workflowQueryTime"
  val CarboniteFreezingTimeMetricPath: NonEmptyList[String] = carboniteMetricsBasePath :+ "freezeTime"

  val CarboniteWorkflowsToCarboniteMetricPath: NonEmptyList[String] = carboniteMetricsBasePath :+ "workflowsToCarbonite"

}
