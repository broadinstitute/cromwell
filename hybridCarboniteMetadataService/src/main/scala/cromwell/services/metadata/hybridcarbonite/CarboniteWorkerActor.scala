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
import scala.util.{Failure, Success, Try}

class CarboniteWorkerActor(carboniterConfig: HybridCarboniteConfig,
                           serviceRegistryActor: ActorRef,
                           ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  val carboniteFreezerActor = context.actorOf(CarbonitingMetadataFreezerActor.props(carboniterConfig, self, serviceRegistryActor, ioActor))

  val backOff = SimpleExponentialBackoff(
    initialInterval = 5.seconds,
    maxInterval = 5.minutes,
    multiplier = 1.1,
    randomizationFactor = 0.0
  )

  scheduleNextCarboniting()


  override def receive: Receive = {
    case s: WorkflowQuerySuccess if s.response.results.nonEmpty => carboniteWorkflow(s.response.results.head.id)
    case _: WorkflowQuerySuccess => scheduleNextCarboniting()
    case f: WorkflowQueryFailure => {
      log.error(f.reason, s"Error while querying workflow to carbonite, will retry.")
      scheduleNextCarboniting()
    }
    case c: CarboniteWorkflowComplete => {
      // after completion of carbonite process reset the backoff
      log.info(s"Carboniting complete for workflow ${c.workflowId}")
      scheduleNextCarboniting(true)
    }
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(carboniteFreezerActor))
    case other => log.error(s"Programmer Error! The CarboniteWorkerActor received unexpected message! ($sender sent $other})")
  }


  def scheduleNextCarboniting(resetBackoff: Boolean = false): Unit = {
    if (resetBackoff) backOff.googleBackoff.reset()
    val duration = Duration(backOff.backoffMillis, MILLISECONDS)
    context.system.scheduler.scheduleOnce(duration)(findWorkflowToCarbonite())
    ()
  }


  def findWorkflowToCarbonite(): Unit = {
    // TODO: [CARBONITE] When carboniting we probably should publish metrics here
    serviceRegistryActor ! QueryForWorkflowsMatchingParameters(CarboniteWorkerActor.findWorkflowToCarboniteQueryParameters)
  }


  def carboniteWorkflow(workflowId: String): Unit = {
    Try(WorkflowId.fromString(workflowId)) match {
      case Success(id: WorkflowId) => {
        log.info(s"Starting carboniting of workflow: $workflowId")
        carboniteFreezerActor ! FreezeMetadata(id)
      }
      case Failure(e) => {
        log.error(e, s"Cannot carbonite workflow $workflowId. Error while converting it to WorkflowId, will retry.")
        scheduleNextCarboniting()
      }
    }
  }
}

object CarboniteWorkerActor {

  sealed trait CarboniteWorkerMessage
  case class CarboniteWorkflowComplete(workflowId: WorkflowId) extends CarboniteWorkerMessage


  val findWorkflowToCarboniteQueryParameters: Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "false",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  )


  def props(carboniterConfig: HybridCarboniteConfig, serviceRegistryActor: ActorRef, ioActor: ActorRef) =
    Props(new CarboniteWorkerActor(carboniterConfig, serviceRegistryActor, ioActor))
}
