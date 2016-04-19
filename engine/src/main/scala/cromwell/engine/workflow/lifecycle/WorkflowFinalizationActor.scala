package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, Props, ActorRef}
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
  case object FinalizationPendingState extends WorkflowFinalizationActorState
  case object FinalizationInProgressState extends WorkflowFinalizationActorState
  case object FinalizationSucceededState extends WorkflowFinalizationActorTerminalState
  case object WorkflowFinalizationFailedState extends WorkflowFinalizationActorTerminalState
  case object FinalizationAbortedState extends WorkflowFinalizationActorTerminalState

  final case class WorkflowFinalizationActorData(backendActors: Map[ActorRef, String],
                                                 successes: Seq[ActorRef],
                                                 failures: Map[ActorRef, Throwable]) extends WorkflowLifecycleActorData

  /**
    * Commands
    */
  sealed trait WorkflowFinalizationActorCommand
  case object StartFinalizationCommand extends WorkflowFinalizationActorCommand
  case object AbortFinalizationCommand extends WorkflowFinalizationActorCommand

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

  override val successState = FinalizationSucceededState
  override val failureState = WorkflowFinalizationFailedState

  override val successResponse = WorkflowFinalizationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowFinalizationFailedResponse(reasons)

  startWith(FinalizationPendingState, WorkflowFinalizationActorData(Map.empty, List.empty, Map.empty))

  when(FinalizationPendingState) {
    case Event(StartFinalizationCommand, _) =>
      val backendFinalizationActors = backendWorkflowActors(workflowDescriptor.backendAssignments)
      if (backendFinalizationActors.isEmpty) {
        sender ! WorkflowFinalizationSucceededResponse // Nothing more to do
        goto(FinalizationSucceededState)
      } else {
        backendFinalizationActors.keys foreach { _ ! Finalize }
        goto(FinalizationInProgressState) using WorkflowFinalizationActorData(backendFinalizationActors, List.empty, Map.empty)
      }
    case Event(AbortFinalizationCommand, _) =>
      context.parent ! WorkflowFinalizationAbortedResponse
      goto(FinalizationAbortedState)
  }

  when(FinalizationInProgressState) {
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
    case Event(AbortFinalizationCommand, _) => ??? // TODO: Handle this
  }

  when(FinalizationSucceededState) { FSM.NullFunction }
  when(WorkflowFinalizationFailedState) { FSM.NullFunction }
  when(FinalizationAbortedState) { FSM.NullFunction }

  /**
    * Makes a BackendWorkflowFinalizationActor for a backend
    */
  override def backendActor(backendName: String, callAssignments: Seq[Call]): Option[ActorRef] = {
    // TODO: Implement this!
    None
  }
}
