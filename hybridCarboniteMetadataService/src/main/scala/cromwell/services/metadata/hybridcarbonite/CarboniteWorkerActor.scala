package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.services.metadata.MetadataArchiveStatus.Unarchived
import cromwell.services.metadata.MetadataService.{QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowComplete
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor.FreezeMetadata
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteWorkerActor(carboniterConfig: HybridCarboniteConfig,
                           serviceRegistryActor: ActorRef,
                           ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  val carboniteFreezerActor = context.actorOf(CarbonitingMetadataFreezerActor.props(carboniterConfig, self, serviceRegistryActor, ioActor))

  val findWorkflowToCarboniteQuery: Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "false",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  )

  val backOff = SimpleExponentialBackoff(
    initialInterval = 5.seconds,
    maxInterval = 5.minutes,
    multiplier = 1.1,
    randomizationFactor = 0.0
  )

  scheduleNextCarboniting(backOff.backoffMillis)


  override def receive: Receive = {
    case s: WorkflowQuerySuccess => {
      if (s.response.results.nonEmpty) carboniteWorkflow(s.response.results.head.id)
      else scheduleNextCarboniting(backOff.backoffMillis)
    }
    case f: WorkflowQueryFailure => {
      log.error(s"Something went wrong while querying to find workflow to carbonite. It will retry in sometime. Error: ${f.reason}")
      scheduleNextCarboniting(backOff.backoffMillis)
    }
    case c: CarboniteWorkflowComplete => {
      // after completion of carbonite process reset the backoff
      log.info(s"Carbonation complete for workflow ${c.workflowId}")
      backOff.googleBackoff.reset()
      scheduleNextCarboniting(backOff.backoffMillis)
    }
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(carboniteFreezerActor))
    case other => log.error(s"Programmer Error! The CarboniteWorkerActor received unexpected message! ($sender sent $other})")
  }


  def scheduleNextCarboniting(interval: Long): Unit = {
    val duration = Duration(interval, MILLISECONDS)
    context.system.scheduler.scheduleOnce(duration)(findWorkflowToCarbonite())
    ()
  }


  def findWorkflowToCarbonite(): Unit = {
    // TODO: [CARBONITE] When carboniting actually does something, we probably metrics here rather than a tick log...
    log.info("Carbonite Worker Tick...")
    serviceRegistryActor ! QueryForWorkflowsMatchingParameters(findWorkflowToCarboniteQuery)
  }


  def carboniteWorkflow(workflowId: String): Unit = {
    log.info(s"Starting Carbonation of workflow: $workflowId")
    carboniteFreezerActor ! FreezeMetadata(WorkflowId.fromString(workflowId))
  }
}

object CarboniteWorkerActor {

  def props(carboniterConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef, ioActor: ActorRef) =
    Props(new CarboniteWorkerActor(carboniterConfig, serviceRegistryActor, ioActor))


  sealed trait CarboniteWorkerMessage
  case class CarboniteWorkflowComplete(workflowId: WorkflowId) extends CarboniteWorkerMessage
}
