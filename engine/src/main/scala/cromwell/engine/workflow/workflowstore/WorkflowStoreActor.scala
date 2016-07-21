package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.services.ServiceRegistryClient
import cromwell.services.metadata.{MetadataValue, MetadataKey, MetadataEvent}
import cromwell.services.metadata.MetadataService.{PutMetadataAction, MetadataPutAcknowledgement}
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success}
import scalaz.NonEmptyList

case class WorkflowStoreActor(store: WorkflowStore) extends LoggingFSM[WorkflowStoreActorState, WorkflowStoreActorData] with ActorLogging with ServiceRegistryClient {

  implicit val ec: ExecutionContext = context.dispatcher

  startWith(Unstarted, WorkflowStoreActorData(None, List.empty))
  self ! InitializerCommand

  when(Unstarted) {
    case Event(InitializerCommand, _) =>
      val work = store.initialize map { unit =>
        log.debug("Workflow store initialization successful")
      }
      addWorkCompletionHooks(InitializerCommand, work)
      goto(Working) using stateData.setCurrentCommand(InitializerCommand, sender)
    case Event(WorkDone, _) =>
      log.error("Shouldn't ever receive {} while Unstarted", WorkDone)
      stay
    case Event(x: WorkflowStoreActorCommand, _) =>
      stay using stateData.withPendingCommand(x, sender)
  }

  when(Idle) {
    case Event(cmd: WorkflowStoreActorCommand, _) =>
      if (stateData.currentOperation.nonEmpty || stateData.pendingOperations.nonEmpty) {
        log.error("Non-empty WorkflowStoreActorData when in Idle state: {}", stateData)
      }
      startNewWork(cmd, sender, stateData.setCurrentCommand(cmd, sender))
  }

  when(Working) {
    case Event(WorkDone, data) =>
      val newData = data.pop
      newData.currentOperation match {
        case None => goto(Idle) using newData
        case Some(WorkflowStoreActorCommandWithSender(cmd, sndr)) => startNewWork(cmd, sndr, newData)
      }
    case Event(cmd: WorkflowStoreActorCommand, data) => stay using data.withPendingCommand(cmd, sender)
  }

  whenUnhandled {
    case Event(MetadataPutAcknowledgement(_), _) =>
      stay // Ignored
    case Event(msg, _) =>
      log.warning("Unexpected message to WorkflowStoreActor in state {} with data {}: {}", stateName, stateData, msg)
      stay
  }

  onTransition {
    case fromState -> toState =>
      log.debug("WorkflowStore moving from {} (using {}) to {} (using {})", fromState, stateData, toState, nextStateData)
  }

  private def startNewWork(command: WorkflowStoreActorCommand, sndr: ActorRef, nextData: WorkflowStoreActorData) = {
    val work: Future[Any] = command match {
      case cmd @ SubmitWorkflow(source) =>
        store.add(NonEmptyList(source)) map { ids =>
          val id = ids.head
          registerIdWithMetadataService(id)
          sndr ! WorkflowSubmittedToStore(id)
          log.info("Workflow {} submitted.", id)
        }
      case cmd @ BatchSubmitWorkflows(sources) =>
        store.add(sources) map { ids =>
          val id = ids.head
          ids foreach registerIdWithMetadataService
          sndr ! WorkflowsBatchSubmittedToStore(ids)
          log.info("Workflows {} submitted.", ids.list.mkString(", "))
        }
      case cmd @ FetchRunnableWorkflows(n) =>
        newWorkflowMessage(n) map { nwm =>
          nwm match {
            case NewWorkflowsToStart(workflows) => log.info("{} new workflows fetched", workflows.size)
            case NoNewWorkflowsToStart => log.debug("No workflows fetched")
            case _ => log.error("Unexpected response from newWorkflowMessage({}): {}", n, nwm)
          }
          sndr ! nwm
        }
      case cmd @ RemoveWorkflow(id) =>
        store.remove(id) map { removed =>
          if (removed) {
            log.debug("Workflow {} removed from store successfully.", id)
          } else {
            log.warning(s"Attempted to remove ID {} from the WorkflowStore but it didn't exist", id)
          }
        }
      case oops =>
        log.error("Unexpected type of start work command: {}", oops.getClass.getSimpleName)
        Future.successful(self ! WorkDone)
    }
    addWorkCompletionHooks(command, work)
    goto(Working) using nextData
  }

  private def addWorkCompletionHooks[A](command: WorkflowStoreActorCommand, work: Future[A]) = {
    work.onComplete {
      case Success(_) =>
        self ! WorkDone
      case Failure(t) =>
        log.error("Error occurred during {}: {} because {}", command.getClass.getSimpleName, t.toString, ExceptionUtils.getStackTrace(t))
        self ! WorkDone
    }
  }

  /**
    * Fetches at most n workflows, and builds the correct response message based on if there were any workflows or not
    */
  private def newWorkflowMessage(maxWorkflows: Int): Future[WorkflowStoreActorResponse] = {
    def fetchRunnableWorkflowsIfNeeded(maxWorkflowsInner: Int, state: StartableState) = {
      if (maxWorkflows > 0) {
        store.fetchRunnableWorkflows(maxWorkflowsInner, state)
      } else {
        Future.successful(List.empty[WorkflowToStart])
      }
    }

    val runnableWorkflows = for {
      restartableWorkflows <- fetchRunnableWorkflowsIfNeeded(maxWorkflows, Restartable)
      submittedWorkflows <- fetchRunnableWorkflowsIfNeeded(maxWorkflows - restartableWorkflows.size, Submitted)
    } yield restartableWorkflows ++ submittedWorkflows

    runnableWorkflows map {
      case x :: xs => NewWorkflowsToStart(NonEmptyList.nel(x, xs))
      case _ => NoNewWorkflowsToStart
    }
  }

  /**
    * Takes the workflow id and sends it over to the metadata service w/ default empty values for inputs/outputs
    */
  private def registerIdWithMetadataService(id: WorkflowId): Unit = {
    val submissionEvents = List(
      MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now.toString)),
      MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
      MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs))
    )

    serviceRegistryActor ! PutMetadataAction(submissionEvents)
  }
}

object WorkflowStoreActor {

  private[workflowstore] case class WorkflowStoreActorCommandWithSender(command: WorkflowStoreActorCommand, sender: ActorRef)

  private[workflowstore] case class WorkflowStoreActorData(currentOperation: Option[WorkflowStoreActorCommandWithSender], pendingOperations: List[WorkflowStoreActorCommandWithSender]) {
    def setCurrentCommand(command: WorkflowStoreActorCommand, sender: ActorRef) = this.copy(currentOperation = Option(WorkflowStoreActorCommandWithSender(command, sender)))
    def withPendingCommand(newCommand: WorkflowStoreActorCommand, sender: ActorRef) = this.copy(pendingOperations = this.pendingOperations :+ WorkflowStoreActorCommandWithSender(newCommand, sender))
    def pop = {
      if (pendingOperations.isEmpty) { WorkflowStoreActorData(None, List.empty) }
      else { WorkflowStoreActorData(Option(pendingOperations.head), pendingOperations.tail) }
    }
  }

  private[workflowstore] sealed trait WorkflowStoreActorState
  private[workflowstore] case object Unstarted extends WorkflowStoreActorState
  private[workflowstore] case object Working extends WorkflowStoreActorState
  private[workflowstore] case object Idle extends WorkflowStoreActorState

  sealed trait WorkflowStoreActorCommand
  case class SubmitWorkflow(source: WorkflowSourceFiles) extends WorkflowStoreActorCommand
  case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFiles]) extends WorkflowStoreActorCommand
  case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorCommand
  case class RemoveWorkflow(id: WorkflowId) extends WorkflowStoreActorCommand

  private case object InitializerCommand extends WorkflowStoreActorCommand
  private case object WorkDone

  sealed trait WorkflowStoreActorResponse
  case class WorkflowSubmittedToStore(workflowId: WorkflowId) extends WorkflowStoreActorResponse
  case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId]) extends WorkflowStoreActorResponse
  case object NoNewWorkflowsToStart extends WorkflowStoreActorResponse
  case class NewWorkflowsToStart(workflows: NonEmptyList[WorkflowToStart]) extends WorkflowStoreActorResponse

  def props(workflowStoreDatabase: WorkflowStore) = Props(WorkflowStoreActor(workflowStoreDatabase)).withDispatcher("akka.dispatchers.api-dispatcher")
}
