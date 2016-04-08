package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor._
import wdl4s.Call

object WorkflowInitializationActor {

  /**
    * States
    */
  sealed trait WorkflowInitializationActorState { def terminal = false }
  sealed trait WorkflowInitializationActorTerminalState extends WorkflowInitializationActorState { override def terminal = true }

  case object InitializationPendingState extends WorkflowInitializationActorState
  case object InitializingWorkflowState extends WorkflowInitializationActorState
  case object InitializationSucceededState extends WorkflowInitializationActorTerminalState
  case object InitializationFailedState extends WorkflowInitializationActorTerminalState
  case object InitializationsAborted extends WorkflowInitializationActorTerminalState

  /**
    * State data:
    *
    * @param initializationActors The mapping from ActorRef to the backend which that actor is initializing.
    * @param successes The list of successful initializations
    * @param failures The map from ActorRef to the reason why that initialization failed.
    */
  final case class WorkflowInitializationActorData(initializationActors: Map[ActorRef, String],
                                                         successes: Seq[ActorRef],
                                                         failures: Map[ActorRef, Throwable])

  /**
    * Commands
    */
  sealed trait WorkflowInitializationActorCommand
  case object StartInitializationCommand extends WorkflowInitializationActorCommand
  case object AbortInitializationCommand extends WorkflowInitializationActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowInitializationActorResponse
  case object WorkflowInitializationSucceededResponse extends WorkflowInitializationActorResponse
  case object WorkflowInitializationAbortedResponse extends WorkflowInitializationActorResponse
  final case class WorkflowInitializationFailedResponse(reasons: Map[ActorRef, Throwable]) extends WorkflowInitializationActorResponse

  def props(workflowTag: String, backendAssignments: Map[Call, String]): Props = Props(new WorkflowInitializationActor(workflowTag, backendAssignments))

}

case class WorkflowInitializationActor(workflowTag: String, backendsToInitialize: Map[Call, String]) extends LoggingFSM[WorkflowInitializationActorState, WorkflowInitializationActorData] {

  startWith(InitializationPendingState, WorkflowInitializationActorData(Map.empty, List.empty, Map.empty))

  when(InitializationPendingState) {
    case Event(StartInitializationCommand, _) =>
      val backendInitializationActors = backendWorkflowInitializationActors // Create each InitializationActor with the list of calls assigned to it.
      backendInitializationActors.keys foreach { _ ! Initialize }
      goto(InitializingWorkflowState) using WorkflowInitializationActorData(backendInitializationActors, List.empty, Map.empty)
    case Event(AbortInitializationCommand, _) =>
      context.parent ! WorkflowInitializationAbortedResponse
      goto(InitializationsAborted)
  }

  when(InitializingWorkflowState) {
    case Event(InitializationSuccess, stateData) =>
      val newData = stateData.copy(successes = stateData.successes ++ List(sender))
      checkForDoneAndTransition(newData)
    case Event(InitializationFailed(reason), stateData) =>
      val newData = stateData.copy(failures = stateData.failures ++ Map(sender -> reason))
      checkForDoneAndTransition(newData)
    case Event(AbortInitializationCommand, _) => ??? // TODO: Handle this
  }

  private def checkForDoneAndTransition(newData: WorkflowInitializationActorData): State = {
    if (checkForDone(newData)) {
      if (newData.failures.isEmpty){
        context.parent ! WorkflowInitializationSucceededResponse
        goto(InitializationSucceededState)
      }
      else {
        context.parent ! WorkflowInitializationFailedResponse(newData.failures)
        goto(InitializationFailedState) using newData
      }
    }
    else {
      stay using newData
    }
  }

  private def checkForDone(stateData: WorkflowInitializationActorData) = {
    // TODO: Check these are actually the same?
    val allActorResponses = stateData.successes ++ stateData.failures.keys
    allActorResponses.size == backendsToInitialize.size
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"Backends initializer for $workflowTag received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> state if state.terminal =>
      log.debug(s"$workflowTag's initializer is terminal. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      log.debug(s"$workflowTag's initializer is transitioning from $fromState to $toState.")
  }

  /**
    * Based on backendAssignments, creates a BackendWorkflowInitializationActor per backend and maps the actor to the
    * name of the backend which it is initializing.
    */
  def backendWorkflowInitializationActors:  Map[ActorRef, String] = ??? //TODO: Implement this!
}
