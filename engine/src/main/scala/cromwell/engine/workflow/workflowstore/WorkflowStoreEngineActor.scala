package cromwell.engine.workflow.workflowstore

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, PoisonPill, Props, Timers}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowProcessingEvents.DescriptionEventValue.PickedUp
import cromwell.core._
import cromwell.core.abort.{WorkflowAbortFailureResponse, WorkflowAbortRequestedResponse, WorkflowAbortedResponse}
import cromwell.engine.instrumentation.WorkflowInstrumentation
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.{WorkflowMetadataHelper, WorkflowProcessingEventPublishing}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowStoreAbortResponse, WorkflowStoreState}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor._
import cromwell.services.instrumentation.CromwellInstrumentationScheduler
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final case class WorkflowStoreEngineActor private(store: WorkflowStore,
                                                  workflowStoreAccess: WorkflowStoreAccess,
                                                  serviceRegistryActor: ActorRef,
                                                  abortAllJobsOnTerminate: Boolean,
                                                  workflowHeartbeatConfig: WorkflowHeartbeatConfig)
  extends LoggingFSM[WorkflowStoreActorState, WorkflowStoreActorData] with ActorLogging with WorkflowInstrumentation with CromwellInstrumentationScheduler with WorkflowMetadataHelper with Timers {

  implicit val ec: ExecutionContext = context.dispatcher

  startWith(Unstarted, WorkflowStoreActorData(None, List.empty))
  self ! InitializerCommand

  val instrumentationAction = () => {
    store.stats map { stats: Map[WorkflowStoreState, Int] =>
      // Update the count for Submitted and Running workflows, defaulting to 0
      val statesMap = stats.withDefault(_ => 0)
      updateWorkflowsQueued(statesMap(WorkflowStoreState.Submitted))
      updateWorkflowsRunning(statesMap(WorkflowStoreState.Running))
      updateWorkflowsOnHold(statesMap(WorkflowStoreState.OnHold))
      updateWorkflowsAborting(statesMap(WorkflowStoreState.Aborting))
    }
    ()
  }

  override def preStart() = {
    startInstrumentationTimer()
    super.preStart()
  }

  override def receive = instrumentationReceive(instrumentationAction).orElse(super.receive)

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
    case Event(ShutdownCommand, _) if abortAllJobsOnTerminate =>
      self ! AbortAllRunningWorkflowsCommandAndStop
      stay()
    case Event(ShutdownCommand, _) =>
      context stop self
      stay()
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
      case FetchRunnableWorkflows(count) =>
        newWorkflowMessage(count) map { response =>
          response match {
            case NewWorkflowsToStart(workflows) =>
              val workflowsIds = workflows.map(_.id).toList
              log.info(
                "{} new workflows fetched by {}: {}",
                workflowsIds.size,
                workflowHeartbeatConfig.cromwellId,
                workflowsIds.mkString(", ")
              )

              workflowsIds foreach { w =>
                WorkflowProcessingEventPublishing.publish(w, workflowHeartbeatConfig.cromwellId, PickedUp, serviceRegistryActor)
              }

            case NoNewWorkflowsToStart => log.debug("No workflows fetched by {}", workflowHeartbeatConfig.cromwellId)
            case _ => log.error("Unexpected response from newWorkflowMessage({}): {}", count, response)
          }
          sndr ! response
        }
      case FindWorkflowsWithAbortRequested(cromwellId) =>
        store.findWorkflowsWithAbortRequested(cromwellId) map {
          ids => sndr ! FindWorkflowsWithAbortRequestedSuccess(ids)
        } recover {
          case t => sndr ! FindWorkflowsWithAbortRequestedFailure(t)
        }
      case AbortWorkflowCommand(id) =>
        store.aborting(id) map { workflowStoreAbortResponse =>
          log.info(s"Abort requested for workflow $id.")
          workflowStoreAbortResponse
        } map {
          case WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted =>
            pushCurrentStateToMetadataService(id, WorkflowAborted)
            sndr ! WorkflowAbortedResponse(id)
          case WorkflowStoreAbortResponse.AbortRequested =>
            pushCurrentStateToMetadataService(id, WorkflowAborting)
            sndr ! WorkflowAbortRequestedResponse(id)
          case WorkflowStoreAbortResponse.NotFound =>
            sndr ! WorkflowAbortFailureResponse(id, new WorkflowNotFoundException(s"Couldn't abort $id because no workflow with that ID is in progress"))
        } recover {
          case t =>
            val message = s"Unable to update workflow store to abort $id"
            log.error(t, message)
            // A generic exception type like RuntimeException will produce a 500 at the API layer, which seems appropriate
            // given we don't know much about what went wrong here.  `t.getMessage` so the cause propagates to the client.
            sndr ! WorkflowAbortFailureResponse(id, new RuntimeException(s"$message: ${t.getMessage}", t))
        }
      case AbortAllRunningWorkflowsCommandAndStop =>
        store.abortAllRunning() map { _ =>
          log.info(s"Aborting all running workflows.")
          self ! PoisonPill
        }
      case WorkflowOnHoldToSubmittedCommand(id) =>
        store.switchOnHoldToSubmitted(id) map { _ =>
          sndr ! WorkflowOnHoldToSubmittedSuccess(id)
          pushCurrentStateToMetadataService(id, WorkflowSubmitted)
          log.info(s"Status changed to 'Submitted' for $id")
        } recover {
          case t =>
            val message = s"Couldn't change the status to 'Submitted' from 'On Hold' for workflow $id"
            log.error(message)
            sndr ! WorkflowOnHoldToSubmittedFailure(id, t)
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
    def fetchStartableWorkflowsIfNeeded = {
      if (maxWorkflows > 0) {
        workflowStoreAccess.fetchStartableWorkflows(maxWorkflows, workflowHeartbeatConfig.cromwellId, workflowHeartbeatConfig.ttl)
      } else {
        Future.successful(List.empty[WorkflowToStart])
      }
    }

    fetchStartableWorkflowsIfNeeded map {
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
  def props(
             workflowStore: WorkflowStore,
             workflowStoreAccess: WorkflowStoreAccess,
             serviceRegistryActor: ActorRef,
             abortAllJobsOnTerminate: Boolean,
             workflowHeartbeatConfig: WorkflowHeartbeatConfig
             ) = {
    Props(WorkflowStoreEngineActor(workflowStore, workflowStoreAccess, serviceRegistryActor, abortAllJobsOnTerminate, workflowHeartbeatConfig)).withDispatcher(EngineDispatcher)
  }

  sealed trait WorkflowStoreEngineActorResponse
  case object NoNewWorkflowsToStart extends WorkflowStoreEngineActorResponse
  final case class NewWorkflowsToStart(workflows: NonEmptyList[WorkflowToStart]) extends WorkflowStoreEngineActorResponse

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

  sealed trait WorkflowOnHoldToSubmittedResponse
  case class WorkflowOnHoldToSubmittedSuccess(workflowId: WorkflowId) extends WorkflowOnHoldToSubmittedResponse
  case class WorkflowOnHoldToSubmittedFailure(workflowId: WorkflowId, failure: Throwable) extends WorkflowOnHoldToSubmittedResponse
}
