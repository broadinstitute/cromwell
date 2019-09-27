package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowSucceeded}
import cromwell.services.metadata.MetadataArchiveStatus.Unarchived
import cromwell.services.metadata.MetadataService.{QueryForWorkflowsMatchingParameters, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowSuccess
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class CarboniteWorkerActor(serviceRegistryActor: ActorRef, ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

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
      log.error(s"Something went wrong while querying to find workflow to carbonite. Error: ${f.reason}")
      scheduleNextCarboniting(backOff.backoffMillis)
    }
    case CarboniteWorkflowSuccess => {
      // it was able to find a workflow to carbonite and after successful carbonation reset the backoff
      backOff.googleBackoff.reset()
      scheduleNextCarboniting(backOff.backoffMillis)
    }
    case ShutdownCommand => stopGracefully()
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


  def carboniteWorkflow(workflowId: String) = {
    log.info(s"Starting Carbonation of workflow: $workflowId")

    // after successful carbonation of a workflow send 'success' command and schedule next carboniting
    self ! CarboniteWorkflowSuccess
  }


  // TODO: [CARBONITE] When the carboniting process is implemented, we might need to implement graceful shutdowns
  def stopGracefully(): Unit = {
    context.stop(self)
  }
}

object CarboniteWorkerActor {

  def props(serviceRegistryActor: ActorRef, ioActor: ActorRef) = Props(new CarboniteWorkerActor(serviceRegistryActor, ioActor))

  case object CarboniteWorkflowSuccess
}
