package cromwell.engine.workflow.lifecycle

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{ActorRef, FSM, OneForOneStrategy, Props}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationSuccess, Finalize}
import cromwell.backend._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.{CallOutputs, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.wdl.WdlTaskCall

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

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, ioActor: ActorRef, jobExecutionMap: JobExecutionMap,
            workflowOutputs: CallOutputs, initializationData: AllBackendInitializationData, copyWorkflowOutputsActor: Option[Props]): Props = {
    Props(new WorkflowFinalizationActor(workflowId, workflowDescriptor, ioActor, jobExecutionMap, workflowOutputs, initializationData, copyWorkflowOutputsActor)).withDispatcher(EngineDispatcher)
  }
}

case class WorkflowFinalizationActor(workflowIdForLogging: WorkflowId,
                                     workflowDescriptor: EngineWorkflowDescriptor,
                                     ioActor: ActorRef,
                                     jobExecutionMap: JobExecutionMap,
                                     workflowOutputs: CallOutputs,
                                     initializationData: AllBackendInitializationData,
                                     copyWorkflowOutputsActorProps: Option[Props])
  extends WorkflowLifecycleActor[WorkflowFinalizationActorState] {

  val tag = self.path.name
  val backendAssignments = workflowDescriptor.backendAssignments

  override val successState = FinalizationSucceededState
  override val failureState = WorkflowFinalizationFailedState

  override def successResponse(data: WorkflowLifecycleActorData) = WorkflowFinalizationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowFinalizationFailedResponse(reasons)

  // If an engine or backend finalization actor (children of this actor) dies, send ourselves the failure and stop the child actor
  override def supervisorStrategy = OneForOneStrategy() {
    case failure => 
      self.tell(FinalizationFailed(failure), sender())
      Stop
  }

  startWith(FinalizationPendingState, WorkflowLifecycleActorData.empty)

  when(FinalizationPendingState) {
    case Event(StartFinalizationCommand, _) =>
      val backendFinalizationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keySet)
          props <- CromwellBackends.backendLifecycleFactoryActorByName(backend).map(
            _.workflowFinalizationActorProps(workflowDescriptor.backendDescriptor, ioActor, calls, filterJobExecutionsForBackend(calls), workflowOutputs, initializationData.get(backend))
          ).get
          actor = context.actorOf(props, backend)
        } yield actor
      }

      val engineFinalizationActor = Try { copyWorkflowOutputsActorProps.map(context.actorOf(_, "CopyWorkflowOutputsActor")).toList }

      val allActors = for {
        backendFinalizationActorsFromTry <- backendFinalizationActors
        engineFinalizationActorFromTry <- engineFinalizationActor
      } yield backendFinalizationActorsFromTry.toList ++ engineFinalizationActorFromTry

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
  private def filterJobExecutionsForBackend(calls: Set[WdlTaskCall]): JobExecutionMap = {
    jobExecutionMap map {
      case (wd, executedKeys) => wd -> (executedKeys filter { jobKey => calls.contains(jobKey.call) })
    } filter {
      case (_, keys) => keys.nonEmpty
    }
  }

  when(FinalizationInProgressState) {
    case Event(FinalizationSuccess, stateData) => checkForDoneAndTransition(stateData.withSuccess(sender))
    case Event(FinalizationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
  }

  when(FinalizationSucceededState) { FSM.NullFunction }
  when(WorkflowFinalizationFailedState) { FSM.NullFunction }
}
