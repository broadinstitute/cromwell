package cromwell.engine.workflow.lifecycle

import akka.actor.{Props, ActorRef}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.Call

import scala.language.postfixOps

object WorkflowFinalizationActor {

  /**
    * States
    */
  sealed trait WorkflowFinalizationActorState extends WorkflowLifecycleActorState
  sealed trait WorkflowFinalizationActorTerminalState extends WorkflowFinalizationActorState
  case object ReadyToFinalize extends WorkflowFinalizationActorState
  case object WorkflowFinalizing extends WorkflowFinalizationActorState
  case object WorkflowFinalizationSucceededState extends WorkflowFinalizationActorTerminalState
  case object WorkflowFinalizationFailedState extends WorkflowFinalizationActorTerminalState
  case object WorkflowFinalizationAbortedState extends WorkflowFinalizationActorTerminalState

  final case class WorkflowFinalizationActorData(backendActors: Map[ActorRef, String],
                                                 successes: Seq[ActorRef],
                                                 failures: Map[ActorRef, Throwable]) extends WorkflowLifecycleActorData

  /**
    * Commands
    */
  sealed trait WorkflowFinalizationActorCommand
  case object StartEngineFinalizationCommand extends WorkflowFinalizationActorCommand
  case object AbortEngineFinalizationCommand extends WorkflowFinalizationActorCommand

  /**
    * Responses
    */
  case object WorkflowFinalizationSucceededResponse extends WorkflowLifecycleSuccessResponse
  case object WorkflowFinalizationAbortedResponse extends WorkflowLifecycleAbortedResponse
  final case class WorkflowFinalizationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(new WorkflowFinalizationActor(workflowId, workflowDescriptor))

}

case class WorkflowFinalizationActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends WorkflowLifecycleActor[WorkflowFinalizationActorState, WorkflowFinalizationActorData] {

  val tag = self.path.name
  val backendAssignments = workflowDescriptor.backendAssignments

  override val successState = WorkflowFinalizationSucceededState
  override val failureState = WorkflowFinalizationFailedState
  override val abortedState = WorkflowFinalizationAbortedState

  override val successResponse = WorkflowFinalizationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowFinalizationFailedResponse(reasons)

  startWith(ReadyToFinalize, WorkflowFinalizationActorData(Map.empty, List.empty, Map.empty))

  when(ReadyToFinalize) {
    case Event(StartEngineFinalizationCommand, _) =>
      // Create each FinalizationActor with the list of calls which ere assigned to that backend.
      val backendFinalizationActors = backendWorkflowActors(backendAssignments)
      backendFinalizationActors.keys foreach { _ ! Finalize }
      goto(WorkflowFinalizing) using WorkflowFinalizationActorData(backendFinalizationActors, List.empty, Map.empty)
    case Event(AbortEngineFinalizationCommand, _) =>
      context.parent ! WorkflowFinalizationAbortedResponse
      goto(WorkflowFinalizationAbortedState)
  }

  when(WorkflowFinalizing) {
    case Event(FinalizationSuccess, stateData) =>
      val newData = stateData.copy(
        backendActors = stateData.backendActors - sender,
        successes = stateData.successes ++ List(sender))
      checkForDoneAndTransition(newData)
    case Event(FinalizationFailed(reason), stateData) =>
      val newData = stateData.copy(
        backendActors = stateData.backendActors - sender,
        failures = stateData.failures ++ Map(sender -> reason))
      checkForDoneAndTransition(newData)
    case Event(AbortEngineFinalizationCommand, stateData) => ??? // TODO: Implement
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.debug(s"$tag state is now terminal. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      log.debug(s"$tag transitioning from $fromState to $toState.")
  }

  /**
    * Makes a BackendWorkflowFinalizationActor for a backend
    */
  override def backendActor(backendName: String, callAssignments: Seq[Call]): Option[ActorRef] = {
    // TODO: Implement this!
    None
  }
}
