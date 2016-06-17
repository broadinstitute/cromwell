package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, Props}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.core.{ExecutionStore, OutputStore, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

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
  final case class WorkflowFinalizationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, executionStore: ExecutionStore,
  outputStore: OutputStore): Props = {
    Props(new WorkflowFinalizationActor(workflowId, workflowDescriptor, executionStore, outputStore))
  }

}

case class WorkflowFinalizationActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor,
                                     executionStore: ExecutionStore, outputStore: OutputStore)
  extends WorkflowLifecycleActor[WorkflowFinalizationActorState] {

  val tag = self.path.name
  val backendAssignments = workflowDescriptor.backendAssignments

  override val successState = FinalizationSucceededState
  override val failureState = WorkflowFinalizationFailedState

  override val successResponse = WorkflowFinalizationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowFinalizationFailedResponse(reasons)

  startWith(FinalizationPendingState, WorkflowLifecycleActorData.empty)

  when(FinalizationPendingState) {
    case Event(StartFinalizationCommand, _) =>
      val backendFinalizationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keys.toSeq)
          props <- CromwellBackends.shadowBackendLifecycleFactory(backend).map(
            _.workflowFinalizationActorProps(workflowDescriptor.backendDescriptor, calls, executionStore, outputStore)
          ).get
          actor = context.actorOf(props)
        } yield actor
      }

      val engineFinalizationActor = Try {
        context.actorOf(CopyWorkflowOutputsActor.props(workflowId, workflowDescriptor, outputStore))
      }

      val allActors = for {
        backendFinalizationActorsFromTry <- backendFinalizationActors
        engineFinalizationActorFromTry <- engineFinalizationActor
      } yield backendFinalizationActorsFromTry.toList.+:(engineFinalizationActorFromTry)

      allActors match {
        case Failure(ex) =>
          sender ! WorkflowFinalizationFailedResponse(Seq(ex))
          goto(WorkflowFinalizationFailedState)
        case Success(actors) if actors.isEmpty =>
          sender ! WorkflowFinalizationSucceededResponse
          goto(FinalizationSucceededState)
        case Success(actors) =>
          val actorSet = actors.toSet
          actorSet.foreach(_ ! Finalize)
          goto(FinalizationInProgressState) using stateData.withActors(actorSet)
        case _ =>
          goto(WorkflowFinalizationFailedState)
      }
  }

  when(FinalizationInProgressState) {
    case Event(FinalizationSuccess, stateData) => checkForDoneAndTransition(stateData.withSuccess(sender))
    case Event(FinalizationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
  }

  when(FinalizationSucceededState) { FSM.NullFunction }
  when(WorkflowFinalizationFailedState) { FSM.NullFunction }
}
