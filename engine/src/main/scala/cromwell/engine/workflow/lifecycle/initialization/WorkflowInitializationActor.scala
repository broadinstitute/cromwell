package cromwell.engine.workflow.lifecycle.initialization

import akka.actor.{ActorRef, FSM, Props}
import common.collections.EnhancedCollections._
import common.exception.AggregatedMessageException
import cromwell.backend.BackendLifecycleActor.BackendActorAbortedResponse
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.backend.{AllBackendInitializationData, BackendWorkflowInitializationActor}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.{PossiblyNotRootWorkflowId, RootWorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import cromwell.engine.workflow.lifecycle.initialization.WorkflowInitializationActor._
import cromwell.engine.workflow.lifecycle.{AbortableWorkflowLifecycleActor, EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}

import scala.util.{Failure, Success, Try}

object WorkflowInitializationActor {

  /**
    * States
    */
  sealed trait WorkflowInitializationActorState extends WorkflowLifecycleActorState
  sealed trait WorkflowInitializationActorTerminalState extends WorkflowInitializationActorState with WorkflowLifecycleActorTerminalState

  case object InitializationPendingState extends WorkflowInitializationActorState
  case object InitializationInProgressState extends WorkflowInitializationActorState
  case object InitializationAbortingState extends WorkflowInitializationActorState
  case object InitializationSucceededState extends WorkflowInitializationActorTerminalState
  case object InitializationFailedState extends WorkflowInitializationActorTerminalState
  case object InitializationsAbortedState extends WorkflowInitializationActorTerminalState

  /**
    * Commands
    */
  sealed trait WorkflowInitializationActorCommand
  case object StartInitializationCommand extends WorkflowInitializationActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowInitializationResponse
  final case class WorkflowInitializationSucceededResponse(initializationData: AllBackendInitializationData) extends WorkflowLifecycleSuccessResponse with WorkflowInitializationResponse
  case object WorkflowInitializationAbortedResponse extends EngineLifecycleActorAbortedResponse with WorkflowInitializationResponse
  final case class WorkflowInitializationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse with WorkflowInitializationResponse

  def props(workflowIdForLogging: PossiblyNotRootWorkflowId,
            rootWorkflowIdForLogging: RootWorkflowId,
            workflowDescriptor: EngineWorkflowDescriptor,
            ioActor: ActorRef,
            serviceRegistryActor: ActorRef,
            restarting: Boolean): Props = {
    Props(new WorkflowInitializationActor(
      workflowIdForLogging = workflowIdForLogging,
      rootWorkflowIdForLogging = rootWorkflowIdForLogging,
      workflowDescriptor = workflowDescriptor,
      ioActor = ioActor,
      serviceRegistryActor = serviceRegistryActor,
      restarting = restarting
    )).withDispatcher(EngineDispatcher)
  }

  case class BackendActorAndBackend(actor: ActorRef, backend: String)
}

case class WorkflowInitializationActor(workflowIdForLogging: PossiblyNotRootWorkflowId,
                                       rootWorkflowIdForLogging: RootWorkflowId,
                                       workflowDescriptor: EngineWorkflowDescriptor,
                                       ioActor: ActorRef,
                                       serviceRegistryActor: ActorRef,
                                       restarting: Boolean)
  extends AbortableWorkflowLifecycleActor[WorkflowInitializationActorState] {

  startWith(InitializationPendingState, WorkflowLifecycleActorData.empty)
  val tag = self.path.name

  override val abortingState = InitializationAbortingState
  override val successState = InitializationSucceededState
  override val failureState = InitializationFailedState
  override val abortedState = InitializationsAbortedState

  override def successResponse(data: WorkflowLifecycleActorData) = {
    val actorsToBackends = backendActorsAndBackends.map(ab => ab.actor -> ab.backend).toMap
    val actorsToData = data.successes.map(ad => ad.actor -> ad.data).toMap
    val allBackendInitializationData = AllBackendInitializationData(actorsToBackends collect { case (a, b) => b -> actorsToData(a) })
    WorkflowInitializationSucceededResponse(allBackendInitializationData)
  }
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowInitializationFailedResponse(reasons)
  override val abortedResponse = WorkflowInitializationAbortedResponse

  private var backendActorsAndBackends: Iterable[BackendActorAndBackend] = _

  when(InitializationPendingState) {
    case Event(StartInitializationCommand, _) =>
      val backendInitializationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).safeMapValues(_.keySet)
          props <- CromwellBackends.backendLifecycleFactoryActorByName(backend).map(factory =>
            factory.workflowInitializationActorProps(workflowDescriptor.backendDescriptor, ioActor, calls, serviceRegistryActor, restarting)
          ).valueOr(errors => throw AggregatedMessageException("Cannot validate backend factories", errors.toList))
          actor = context.actorOf(props, backend)
        } yield BackendActorAndBackend(actor, backend)
      }

      backendInitializationActors match {
        case Failure(ex) =>
          sender() ! WorkflowInitializationFailedResponse(Seq(ex))
          goto(InitializationFailedState)
        case Success(actors) if actors.isEmpty =>
          backendActorsAndBackends = List.empty
          sender() ! WorkflowInitializationSucceededResponse(AllBackendInitializationData.empty)
          goto(InitializationSucceededState)
        case Success(actors) =>
          backendActorsAndBackends = actors
          val actorSet: Set[ActorRef] = actors.map(_.actor).toSet
          actorSet.foreach(_ ! Initialize)
          goto(InitializationInProgressState) using stateData.withActors(actorSet)
      }

    case Event(EngineLifecycleActorAbortCommand, _) =>
      context.parent ! WorkflowInitializationAbortedResponse
      goto(InitializationsAbortedState)
  }

  when(InitializationInProgressState) {
    case Event(InitializationSuccess(initData), stateData) => checkForDoneAndTransition(stateData.withSuccess(sender(), initData))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender(), reason))
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      stateData.actors foreach { _ ! BackendWorkflowInitializationActor.Abort }
      goto(InitializationAbortingState)
  }

  when(InitializationAbortingState) {
    case Event(InitializationSuccess(initData), stateData) => checkForDoneAndTransition(stateData.withSuccess(sender(), initData))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender(), reason))
    case Event(BackendActorAbortedResponse, stateData) => checkForDoneAndTransition(stateData.withAborted(sender()))
  }

  when(InitializationSucceededState) { FSM.NullFunction }
  when(InitializationFailedState) { FSM.NullFunction }
  when(InitializationsAbortedState) { FSM.NullFunction }
}
