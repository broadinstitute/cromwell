package cromwell.engine.workflow.workflowstore

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowId
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.{WorkflowStoreActorState, _}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final case class WorkflowStoreEngineActor(store: WorkflowStore, serviceRegistryActor: ActorRef)
  extends LoggingFSM[WorkflowStoreActorState, WorkflowStoreActorData] with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  startWith(Unstarted, WorkflowStoreActorData(None, List.empty))
  self ! InitializerCommand

  when(Unstarted) {
    case Event(InitializerCommand, _) =>
      val work = store.initialize map { _ =>
        log.debug("Workflow store initialization successful")
      }
      addWorkCompletionHooks(InitializerCommand, work)
      goto(Working) using stateData.withCurrentCommand(InitializerCommand, sender)
    case Event(x: WorkflowStoreActorEngineCommand, _) =>
      stay using stateData.withPendingCommand(x, sender)
  }

  when(Idle) {
    case Event(cmd: WorkflowStoreActorEngineCommand, _) =>
      if (stateData.currentOperation.nonEmpty || stateData.pendingOperations.nonEmpty) {
        log.error("Non-empty WorkflowStoreActorData when in Idle state: {}", stateData)
      }
      startNewWork(cmd, sender, stateData.withCurrentCommand(cmd, sender))
  }

  when(Working) {
    case Event(WorkDone, data) =>
      val newData = data.pop
      newData.currentOperation match {
        case None => goto(Idle) using newData
        case Some(WorkflowStoreActorCommandWithSender(cmd, sndr)) => startNewWork(cmd, sndr, newData)
      }
    case Event(cmd: WorkflowStoreActorEngineCommand, data) => stay using data.withPendingCommand(cmd, sender)
  }

  whenUnhandled {
    case Event(msg, _) =>
      log.warning("Unexpected message to WorkflowStoreActor in state {} with data {}: {}", stateName, stateData, msg)
      stay
  }

  onTransition {
    case fromState -> toState =>
      log.debug("WorkflowStore moving from {} (using {}) to {} (using {})", fromState, stateData, toState, nextStateData)
  }

  private def startNewWork(command: WorkflowStoreActorEngineCommand, sndr: ActorRef, nextData: WorkflowStoreActorData) = {
    val work: Future[Any] = command match {
      case cmd @ FetchRunnableWorkflows(n) =>
        newWorkflowMessage(n) map { nwm =>
          nwm match {
            case NewWorkflowsToStart(workflows) => log.info("{} new workflows fetched", workflows.toList.size)
            case NoNewWorkflowsToStart => log.debug("No workflows fetched")
            case _ => log.error("Unexpected response from newWorkflowMessage({}): {}", n, nwm)
          }
          sndr ! nwm
        }
      case cmd @ AbortWorkflow(id, manager) =>
        store.remove(id) map { removed =>
          if (removed) {
            manager ! WorkflowManagerActor.AbortWorkflowCommand(id, sndr)
            log.debug(s"Workflow $id removed from the workflow store, abort requested.")
          } else {
            sndr ! WorkflowAbortFailed(id, new WorkflowNotFoundException(s"Couldn't abort $id because no workflow with that ID is in progress"))
          }
        } recover {
          case t =>
            val message = s"Error aborting workflow $id: could not remove from workflow store"
            log.error(t, message)
            // A generic exception type like RuntimeException will produce a 500 at the API layer, which seems appropriate
            // given we don't know much about what went wrong here.  `t.getMessage` so the cause propagates to the client.
            val e = new RuntimeException(s"$message: ${t.getMessage}", t)
            sndr ! WorkflowAbortFailed(id, e)
        }
      case cmd @ RemoveWorkflow(id) =>
        store.remove(id) map { removed =>
          if (removed) {
            log.debug("Workflow {} removed from store successfully.", id)
          } else {
            log.warning(s"Attempted to remove ID {} from the WorkflowStore but it didn't exist", id)
          }
        } recover {
          case t =>
            log.error(t, s"Unable to remove workflow $id from workflow store")
        }
      case oops =>
        log.error("Unexpected type of start work command: {}", oops.getClass.getSimpleName)
        Future.successful(self ! WorkDone)
    }
    addWorkCompletionHooks(command, work)
    goto(Working) using nextData
  }

  private def addWorkCompletionHooks[A](command: WorkflowStoreActorEngineCommand, work: Future[A]) = {
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
  private def newWorkflowMessage(maxWorkflows: Int): Future[WorkflowStoreEngineActorResponse] = {
    def fetchRunnableWorkflowsIfNeeded(maxWorkflowsInner: Int, state: StartableState) = {
      if (maxWorkflows > 0) {
        store.fetchRunnableWorkflows(maxWorkflowsInner, state)
      } else {
        Future.successful(List.empty[WorkflowToStart])
      }
    }

    val runnableWorkflows = for {
      restartableWorkflows <- fetchRunnableWorkflowsIfNeeded(maxWorkflows, WorkflowStoreState.Restartable)
      submittedWorkflows <- fetchRunnableWorkflowsIfNeeded(maxWorkflows - restartableWorkflows.size, WorkflowStoreState.Submitted)
    } yield restartableWorkflows ++ submittedWorkflows

    runnableWorkflows map {
      case x :: xs => NewWorkflowsToStart(NonEmptyList.of(x, xs: _*))
      case _ => NoNewWorkflowsToStart
    } recover {
      case e =>
        // Log the error but return a successful Future so as not to hang future workflow store polls.
        log.error(e, "Error trying to fetch new workflows")
        NoNewWorkflowsToStart
    }
  }
}

object WorkflowStoreEngineActor {
  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreEngineActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(EngineDispatcher)
  }

  sealed trait WorkflowStoreEngineActorResponse
  case object NoNewWorkflowsToStart extends WorkflowStoreEngineActorResponse
  final case class NewWorkflowsToStart(workflows: NonEmptyList[WorkflowToStart]) extends WorkflowStoreEngineActorResponse
  final case class WorkflowAborted(workflowId: WorkflowId) extends WorkflowStoreEngineActorResponse
  final case class WorkflowAbortFailed(workflowId: WorkflowId, reason: Throwable) extends WorkflowStoreEngineActorResponse


  final case class WorkflowStoreActorCommandWithSender(command: WorkflowStoreActorEngineCommand, sender: ActorRef)

  final case class WorkflowStoreActorData(currentOperation: Option[WorkflowStoreActorCommandWithSender], pendingOperations: List[WorkflowStoreActorCommandWithSender]) {
    def withCurrentCommand(command: WorkflowStoreActorEngineCommand, sender: ActorRef) = this.copy(currentOperation = Option(WorkflowStoreActorCommandWithSender(command, sender)))
    def withPendingCommand(newCommand: WorkflowStoreActorEngineCommand, sender: ActorRef) = this.copy(pendingOperations = this.pendingOperations :+ WorkflowStoreActorCommandWithSender(newCommand, sender))
    def pop = {
      if (pendingOperations.isEmpty) { WorkflowStoreActorData(None, List.empty) }
      else { WorkflowStoreActorData(Option(pendingOperations.head), pendingOperations.tail) }
    }
  }

  sealed trait WorkflowStoreActorState
  case object Unstarted extends WorkflowStoreActorState
  case object Working extends WorkflowStoreActorState
  case object Idle extends WorkflowStoreActorState
}
