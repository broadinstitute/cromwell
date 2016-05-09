package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, ActorRef, Props}
import cromwell.backend.BackendLifecycleActor.BackendActorAbortedResponse
import cromwell.backend.BackendWorkflowInitializationActor
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.BackendWorkflowInitializationActor._
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor._
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.Call

import scala.language.postfixOps
import scala.util.{Success, Failure, Try}

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
  case object WorkflowInitializationSucceededResponse extends WorkflowLifecycleSuccessResponse
  case object WorkflowInitializationAbortedResponse extends EngineLifecycleActorAbortedResponse
  final case class WorkflowInitializationFailedResponse(reasons: Seq[Throwable]) extends WorkflowLifecycleFailureResponse

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(new WorkflowInitializationActor(workflowId, workflowDescriptor))
}

case class WorkflowInitializationActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends WorkflowLifecycleActor[WorkflowInitializationActorState] {

  startWith(InitializationPendingState, WorkflowLifecycleActorData.empty)
  val tag = self.path.name

  override val abortingState = InitializationAbortingState
  override val successState = InitializationSucceededState
  override val failureState = InitializationFailedState
  override val abortedState = InitializationsAbortedState

  override val successResponse = WorkflowInitializationSucceededResponse
  override def failureResponse(reasons: Seq[Throwable]) = WorkflowInitializationFailedResponse(reasons)
  override val abortedResponse = WorkflowInitializationAbortedResponse

  when(InitializationPendingState) {
    case Event(StartInitializationCommand, _) =>
      val backendInitializationActors = Try {
        for {
          (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keys.toSeq)
          backendConfiguration = BackendConfiguration.backendConfigurationDescriptor(backend).get
          props <- CromwellBackends.shadowBackendLifecycleFactory(backend).map(factory =>
            factory.workflowInitializationActorProps(workflowDescriptor.backendDescriptor, calls, backendConfiguration)
          ).get
          actor = context.actorOf(props)
        } yield (actor, backend)
      }

      backendInitializationActors match {
        case Failure(ex) =>
          sender ! WorkflowInitializationFailedResponse(Seq(ex))
          goto(InitializationFailedState)
        case Success(actors) if actors.isEmpty=>
          sender ! WorkflowInitializationSucceededResponse
          goto(InitializationSucceededState)
        case Success(actors) =>
          actors.keys.foreach(_ ! Initialize)
          goto(InitializationInProgressState) using stateData.withBackendActors(actors)
      }

    case Event(InitializationAbortingState, _) =>
      context.parent ! WorkflowInitializationAbortedResponse
      goto(InitializationsAbortedState)
  }

  when(InitializationInProgressState) {
    case Event(InitializationSuccess, stateData) => checkForDoneAndTransition(stateData.withSuccess(sender))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      stateData.backendActors.keys foreach { _ ! BackendWorkflowInitializationActor.Abort }
      goto(InitializationAbortingState)
  }

  when(InitializationAbortingState) {
    case Event(InitializationSuccess, stateData) => checkForDoneAndTransition(stateData.withSuccess(sender))
    case Event(InitializationFailed(reason), stateData) => checkForDoneAndTransition(stateData.withFailure(sender, reason))
    case Event(BackendActorAbortedResponse, stateData) => checkForDoneAndTransition(stateData.withAborted(sender))
  }

  when(InitializationSucceededState) { FSM.NullFunction }
  when(InitializationFailedState) { FSM.NullFunction }
  when(InitializationsAbortedState) { FSM.NullFunction }
}
