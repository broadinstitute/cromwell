package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.Call

import scala.language.postfixOps

object WorkflowInitializationActor {

  /**
    * States
    */
  sealed trait WorkflowInitializationActorState extends WorkflowLifecycleActorState
  sealed trait WorkflowInitializationActorTerminalState extends WorkflowInitializationActorState with WorkflowLifecycleActorTerminalState

  case object InitializationPendingState extends WorkflowInitializationActorState
  case object InitializingWorkflowState extends WorkflowInitializationActorState
  case object InitializationSucceededState extends WorkflowInitializationActorTerminalState
  case object InitializationFailedState extends WorkflowInitializationActorTerminalState
  case object InitializationsAbortedState extends WorkflowInitializationActorTerminalState

  /**
    * State data:
    *
    * @param backendActors The mapping from ActorRef to the backend which that actor is initializing.
    * @param successes The list of successful initializations
    * @param failures The map from ActorRef to the reason why that initialization failed.
    */
  final case class WorkflowInitializationActorData(backendActors: Map[ActorRef, String],
                                                   successes: Seq[ActorRef],
                                                   failures: Map[ActorRef, Throwable]) extends WorkflowLifecycleActorData

  /**
    * Commands
    */
  sealed trait WorkflowInitializationActorCommand
  case object StartInitializationCommand extends WorkflowInitializationActorCommand
  case object AbortInitializationCommand extends WorkflowInitializationActorCommand

  /**
    * Responses
    */
  case object WorkflowInitializationSucceededResponse extends WorkflowLifecycleSuccessResponse
  case object WorkflowInitializationAbortedResponse extends WorkflowLifecycleAbortedResponse
  final case class WorkflowInitializationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(new WorkflowInitializationActor(workflowId, workflowDescriptor))
}

case class WorkflowInitializationActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends WorkflowLifecycleActor[WorkflowInitializationActorState, WorkflowInitializationActorData] {

  startWith(InitializationPendingState, WorkflowInitializationActorData(Map.empty, List.empty, Map.empty))
  val tag = self.path.name

  override val successState = InitializationSucceededState
  override val failureState = InitializationFailedState
  override val abortedState = InitializationFailedState

  override val successResponse = WorkflowInitializationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowInitializationFailedResponse(reasons)

  when(InitializationPendingState) {
    case Event(StartInitializationCommand, _) =>
      val backendInitializationActors = backendWorkflowActors(workflowDescriptor.backendAssignments)
      if (backendInitializationActors.isEmpty) {
        sender ! WorkflowInitializationSucceededResponse // Nothing more to do
        goto(InitializationFailedState)
      } else {
        backendInitializationActors.keys foreach { _ ! Initialize }
        goto(InitializingWorkflowState) using WorkflowInitializationActorData(backendInitializationActors, List.empty, Map.empty)
      }
    case Event(AbortInitializationCommand, _) =>
      context.parent ! WorkflowInitializationAbortedResponse
      goto(InitializationsAbortedState)
  }

  when(InitializingWorkflowState) {
    case Event(InitializationSuccess, stateData) =>
      val newData = stateData.copy(
        backendActors = stateData.backendActors - sender,
        successes = stateData.successes ++ List(sender))
      checkForDoneAndTransition(newData)
    case Event(InitializationFailed(reason), stateData) =>
      val newData = stateData.copy(
        backendActors = stateData.backendActors - sender,
        failures = stateData.failures ++ Map(sender -> reason))
      checkForDoneAndTransition(newData)
    case Event(AbortInitializationCommand, _) => ??? // TODO: Handle this
  }

  /**
    * Makes a BackendWorkflowInitializationActor for a backend
    */
  override def backendActor(backendName: String, callAssignments: Seq[Call]): Option[ActorRef] = {
    // TODO: Implement this!
    None
  }
}
