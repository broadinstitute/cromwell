package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, FSM, Props}
import cromwell.backend.BackendLifecycleActor.BackendActorAbortedResponse
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.backend.{AllBackendInitializationData, BackendWorkflowInitializationActor}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._

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
  final case class WorkflowInitializationSucceededResponse(initializationData: AllBackendInitializationData) extends WorkflowLifecycleSuccessResponse
  case object WorkflowInitializationAbortedResponse extends EngineLifecycleActorAbortedResponse
  final case class WorkflowInitializationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse

  def props(workflowId: WorkflowId,
            workflowDescriptor: EngineWorkflowDescriptor,
            serviceRegistryActor: ActorRef): Props = {
    Props(new WorkflowInitializationActor(workflowId, workflowDescriptor, serviceRegistryActor)).withDispatcher(EngineDispatcher)
  }

  case class BackendActorAndBackend(actor: ActorRef, backend: String)
}

case class WorkflowInitializationActor(workflowIdForLogging: WorkflowId,
                                       workflowDescriptor: EngineWorkflowDescriptor,
                                       serviceRegistryActor: ActorRef)
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
    val allBackendInitializationData = AllBackendInitializationData(actorsToBackends collect { case (a, b) => b -> actorsToData.get(a).get })
    WorkflowInitializationSucceededResponse(allBackendInitializationData)
  }
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowInitializationFailedResponse(reasons)
  override val abortedResponse = WorkflowInitializationAbortedResponse

  private var backendActorsAndBackends: Traversable[BackendActorAndBackend] = _

  when(InitializationPendingState) {
    case Event(StartInitializationCommand, _) =>
      val backendInitializationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keySet)
          props <- CromwellBackends.backendLifecycleFactoryActorByName(backend).map(factory =>
            factory.workflowInitializationActorProps(workflowDescriptor.backendDescriptor, calls, serviceRegistryActor)
          ).get
          actor = context.actorOf(props, backend)
        } yield BackendActorAndBackend(actor, backend)
      }

      backendInitializationActors match {
        case Failure(ex) =>
          sender ! WorkflowInitializationFailedResponse(Seq(ex))
          goto(InitializationFailedState)
        case Success(actors) if actors.isEmpty =>
          backendActorsAndBackends = List.empty
          sender ! WorkflowInitializationSucceededResponse(AllBackendInitializationData.empty)
          goto(InitializationSucceededState)
        case Success(actors) =>
          backendActorsAndBackends = actors
          val actorSet: Set[ActorRef] = actors.map(_.actor).toSet
          actorSet.foreach(_ ! Initialize)
          goto(InitializationInProgressState) using stateData.withActors(actorSet)
      }

    case Event(InitializationAbortingState, _) =>
      context.parent ! WorkflowInitializationAbortedResponse
      goto(InitializationsAbortedState)
  }

  when(InitializationInProgressState) {
    case Event(InitializationSuccess(initData), stateData) => checkForDoneAndTransition(stateData.withSuccess(sender, initData))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      stateData.actors foreach { _ ! BackendWorkflowInitializationActor.Abort }
      goto(InitializationAbortingState)
  }

  when(InitializationAbortingState) {
    case Event(InitializationSuccess(initData), stateData) => checkForDoneAndTransition(stateData.withSuccess(sender, initData))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
    case Event(BackendActorAbortedResponse, stateData) => checkForDoneAndTransition(stateData.withAborted(sender))
  }

  when(InitializationSucceededState) { FSM.NullFunction }
  when(InitializationFailedState) { FSM.NullFunction }
  when(InitializationsAbortedState) { FSM.NullFunction }
}
