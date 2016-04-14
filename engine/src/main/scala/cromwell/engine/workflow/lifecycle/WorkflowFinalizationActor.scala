package cromwell.engine.workflow.lifecycle

import akka.actor.{Props, ActorRef, LoggingFSM}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor._
import wdl4s.Call

object WorkflowFinalizationActor {

  /**
    * States
    */
  sealed trait WorkflowFinalizationActorState { def terminal = false }
  sealed trait WorkflowFinalizationActorTerminalState extends WorkflowFinalizationActorState { override def terminal = true }
  case object ReadyToFinalize extends WorkflowFinalizationActorState
  case object EngineFinalizing extends WorkflowFinalizationActorState
  case object EngineFinalizationSucceededState extends WorkflowFinalizationActorTerminalState
  case object EngineFinalizationFailedState extends WorkflowFinalizationActorTerminalState
  case object EngineFinalizationAbortedState extends WorkflowFinalizationActorTerminalState

  /**
    * State data:
    *
    * @param finalizationActors The mapping from ActorRef to the name of the backend which that actor is finalizing.
    * @param successes The list of successful initializations
    * @param failures The map from ActorRef to the reason why that initialization failed.
    */
  final case class WorkflowFinalizationActorData(finalizationActors: Map[ActorRef, String],
                                                       successes: Seq[ActorRef],
                                                       failures: Map[ActorRef, Throwable])

  /**
    * Commands
    */
  sealed trait WorkflowFinalizationActorCommand
  case object StartEngineFinalizationCommand extends WorkflowFinalizationActorCommand
  case object AbortEngineFinalizationCommand extends WorkflowFinalizationActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowFinalizationActorResponse
  case object WorkflowFinalizationSucceededResponse extends WorkflowFinalizationActorResponse
  case object WorkflowFinalizationAbortedResponse extends WorkflowFinalizationActorResponse
  final case class WorkflowFinalizationFailedResponse(reasons: Seq[Throwable]) extends WorkflowFinalizationActorResponse

  def props(workflowTag: String, backendAssignments: Map[Call, String]): Props = Props(new WorkflowFinalizationActor(workflowTag, backendAssignments))

}

case class WorkflowFinalizationActor(workflowTag: String, backendsToFinalize: Map[Call, String]) extends LoggingFSM[WorkflowFinalizationActorState, WorkflowFinalizationActorData] {

  val tag = s"$workflowTag finalizer"
  startWith(ReadyToFinalize, WorkflowFinalizationActorData(Map.empty, List.empty, Map.empty))

  when(ReadyToFinalize) {
    case Event(StartEngineFinalizationCommand, _) =>
      // Create each FinalizationActor with the list of calls which ere assigned to that backend.
      val backendInitializationActors = backendWorkflowFinalizationActors
      backendInitializationActors.keys foreach { _ ! Finalize }
      goto(EngineFinalizing) using WorkflowFinalizationActorData(backendInitializationActors, List.empty, Map.empty)
    case Event(AbortEngineFinalizationCommand, _) =>
      context.parent ! WorkflowFinalizationAbortedResponse
      goto(EngineFinalizationAbortedState)
  }

  when(EngineFinalizing) {
    case Event(FinalizationSuccess, stateData) =>
      val newData = stateData.copy(successes = stateData.successes ++ List(sender))
      checkForDoneAndTransition(newData)
    case Event(FinalizationFailed(reason), stateData) =>
      val newData = stateData.copy(failures = stateData.failures ++ Map(sender -> reason))
      checkForDoneAndTransition(newData)
    case Event(AbortEngineFinalizationCommand, stateData) => ??? // TODO: Implement
  }

  private def checkForDoneAndTransition(newData: WorkflowFinalizationActorData): State = {
    if (checkForDone(newData)) {
      if (newData.failures.isEmpty) {
        context.parent ! WorkflowFinalizationSucceededResponse
        goto(EngineFinalizationSucceededState) using newData
      }
      else {
        context.parent ! WorkflowFinalizationFailedResponse(newData.failures.values.toSeq)
        goto(EngineFinalizationFailedState)
      }
    }
    else {
      stay using newData
    }
  }

  private def checkForDone(stateData: WorkflowFinalizationActorData) = {
    // TODO: Check these are actually the same?
    val allActorResponses = stateData.successes ++ stateData.failures.keys
    allActorResponses.size == backendsToFinalize.size
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"Backends finalizer for $workflowTag received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.debug(s"$tag is done. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      log.debug(s"$tag is transitioning from $fromState to $toState.")
  }

  /**
    * Based on the backendAssignments, creates a BackendWorkflowFinalizationActor per backend and maps the actor to the
    * name of the backend which it is finalizing.
    */
  def backendWorkflowFinalizationActors:  Map[ActorRef, String] = ??? //TODO: Implement this!
}
