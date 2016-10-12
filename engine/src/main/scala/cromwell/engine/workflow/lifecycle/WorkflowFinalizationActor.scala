package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, Props}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.backend._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.{CallOutputs, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.TaskCall

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

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap,
            workflowOutputs: CallOutputs, initializationData: AllBackendInitializationData): Props = {
    Props(new WorkflowFinalizationActor(workflowId, workflowDescriptor, jobExecutionMap, workflowOutputs, initializationData)).withDispatcher(EngineDispatcher)
  }
}

case class WorkflowFinalizationActor(workflowIdForLogging: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor,
                                     jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs, initializationData: AllBackendInitializationData)
  extends WorkflowLifecycleActor[WorkflowFinalizationActorState] {

  val tag = self.path.name
  val backendAssignments = workflowDescriptor.backendAssignments

  override val successState = FinalizationSucceededState
  override val failureState = WorkflowFinalizationFailedState

  override def successResponse(data: WorkflowLifecycleActorData) = WorkflowFinalizationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowFinalizationFailedResponse(reasons)

  startWith(FinalizationPendingState, WorkflowLifecycleActorData.empty)

  when(FinalizationPendingState) {
    case Event(StartFinalizationCommand, _) =>
      val backendFinalizationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keySet)
          props <- CromwellBackends.backendLifecycleFactoryActorByName(backend).map(
            _.workflowFinalizationActorProps(workflowDescriptor.backendDescriptor, calls, filterJobExecutionsForBackend(calls), workflowOutputs, initializationData.get(backend))
          ).get
          actor = context.actorOf(props, backend)
        } yield actor
      }

      val engineFinalizationActor = Try {
        context.actorOf(CopyWorkflowOutputsActor.props(workflowIdForLogging, workflowDescriptor, workflowOutputs, initializationData),
          "CopyWorkflowOutputsActor")
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
  
  // Only send to each backend the jobs that it executed
  private def filterJobExecutionsForBackend(calls: Set[TaskCall]): JobExecutionMap = {
    jobExecutionMap map {
      case (wd, executedKeys) => wd -> (executedKeys filter { jobKey => calls.contains(jobKey.call) })
    } filter {
      case (wd, keys) => keys.nonEmpty
    }
  }

  when(FinalizationInProgressState) {
    case Event(FinalizationSuccess, stateData) => checkForDoneAndTransition(stateData.withSuccess(sender))
    case Event(FinalizationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
  }

  when(FinalizationSucceededState) { FSM.NullFunction }
  when(WorkflowFinalizationFailedState) { FSM.NullFunction }
}
