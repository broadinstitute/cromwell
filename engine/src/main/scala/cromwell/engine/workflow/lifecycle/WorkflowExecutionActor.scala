package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.engine.backend.BackendCallJobDescriptor
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor._
import wdl4s.Call

object WorkflowExecutionActor {

  /**
    * States
    */
  sealed trait WorkflowExecutionActorState { def terminal = false }
  sealed trait WorkflowExecutionActorTerminalState extends WorkflowExecutionActorState { override def terminal = true }

  case object WorkflowExecutionPendingState extends WorkflowExecutionActorState
  case object WorkflowExecutionInProgressState extends WorkflowExecutionActorState
  case object WorkflowExecutionSuccessfulState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionFailedState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionAbortedState extends WorkflowExecutionActorTerminalState

  /**
    * State data
    */
  final case class WorkflowExecutionActorData()

  /**
    * Commands
    */
  sealed trait WorkflowExecutionActorCommand
  case object StartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object RestartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object AbortExecutingWorkflowCommand extends WorkflowExecutionActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowExecutionActorResponse
  case object WorkflowExecutionSucceededResponse extends WorkflowExecutionActorResponse
  case object WorkflowExecutionAbortedResponse extends WorkflowExecutionActorResponse
  final case class WorkflowExecutionFailedResponse(reasons: Map[ActorRef, Throwable]) extends WorkflowExecutionActorResponse

  def props(workflowTag: String, backendAssignments: Map[Call, String]): Props = Props(WorkflowExecutionActor(workflowTag, backendAssignments))
}

final case class WorkflowExecutionActor(workflowTag: String, backendAssignments: Map[Call, String]) extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] {

  val tag = s"$workflowTag executor"
  startWith(WorkflowExecutionPendingState, WorkflowExecutionActorData())

  when(WorkflowExecutionPendingState) {
    case Event(StartExecutingWorkflowCommand, _) =>
      // TODO: Start executing
      goto(WorkflowExecutionInProgressState)
    case Event(RestartExecutingWorkflowCommand, _) =>
      // TODO: Restart executing
      goto(WorkflowExecutionInProgressState)
    case Event(AbortExecutingWorkflowCommand, _) =>
      context.parent ! WorkflowExecutionAbortedResponse
      goto(WorkflowExecutionAbortedState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(AbortExecutingWorkflowCommand, stateData) => ??? // TODO: Implement!
    case Event(_, _) => ??? // TODO: Lots of extra stuff to include here...
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.debug(s"$tag done. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      log.debug(s"$tag transitioning from $fromState to $toState.")
  }

  /**
    * Creates an appropraite BackendJobExecutionActor to run a single job, according to the backend assignments.
    */
  def backendForExecution(jobDescriptor: BackendCallJobDescriptor): ActorRef = ??? //TODO: Implement!
}
