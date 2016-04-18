package cromwell.engine.workflow

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.{EngineWorkflowDescriptor, WorkflowSourceFiles}

import cromwell.engine.workflow.lifecycle.ShadowMaterializeWorkflowDescriptorActor.{ShadowMaterializeWorkflowDescriptorFailureResponse, ShadowMaterializeWorkflowDescriptorSuccessResponse, ShadowMaterializeWorkflowDescriptorCommand}
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.{WorkflowExecutionSucceededResponse, WorkflowExecutionFailedResponse, RestartExecutingWorkflowCommand, StartExecutingWorkflowCommand}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{WorkflowFinalizationSucceededResponse, WorkflowFinalizationFailedResponse, StartEngineFinalizationCommand}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse, StartInitializationCommand}
import cromwell.engine.workflow.lifecycle.{ShadowMaterializeWorkflowDescriptorActor, WorkflowFinalizationActor, WorkflowExecutionActor, WorkflowInitializationActor}
import cromwell.engine.workflow.ShadowWorkflowActor._
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowOutputs

object ShadowWorkflowActor {

  /**
    * Commands to the ShadowWorkflowActor
    */
  sealed trait ShadowWorkflowActorCommand
  case object StartWorkflowCommand extends ShadowWorkflowActorCommand
  case object AbortWorkflowCommand extends ShadowWorkflowActorCommand

  /**
    * Responses from the ShadowWorkflowActor
    */
  sealed trait ShadowWorkflowActorResponse
  case class ShadowWorkflowSucceededResponse(workflowId: WorkflowId, outputs: WorkflowOutputs) extends ShadowWorkflowActorResponse
  case class ShadowWorkflowFailedResponse(workflowId: WorkflowId, inState: ShadowWorkflowActorState, reasons: Seq[Throwable]) extends Exception with ShadowWorkflowActorResponse

  /**
    * States for the ShadowWorkflowActor FSM
    */
  sealed trait ShadowWorkflowActorState { def terminal = false }
  sealed trait ShadowWorkflowActorTerminalState extends ShadowWorkflowActorState { override def terminal = true}

  /**
    * Waiting for a Start or Restart command.
    */
  case object WorkflowUnstartedState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor is created. It now needs to validate and create its WorkflowDescriptor
    */
  case object MaterializingWorkflowDescriptorState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has a descriptor. It now needs to create its set of WorkflowBackendActors
    */
  case object InitializingWorkflowState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has a descriptor and a set of backends. It now needs to execute the hell out of its jobs.
    */
  case object ExecutingWorkflowState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has completed. So we're now finalizing whatever needs to be finalized on the backends
    */
  case object FinalizingWorkflowState extends ShadowWorkflowActorState

  /**
    * We're done.
    */
  case object WorkflowSucceededState extends ShadowWorkflowActorTerminalState

  /**
    * We're done. But we were aborted in the middle of the lifecycle.
    */
  case object WorkflowAbortedState extends ShadowWorkflowActorTerminalState

  /**
    * We're done. But the workflow has failed.
    */
  case object WorkflowFailedState extends ShadowWorkflowActorTerminalState

  sealed trait ShadowWorkflowActorData
  case object EmptyShadowWorkflowActorData extends ShadowWorkflowActorData
  case class ShadowWorkflowActorDataWithDescriptor(workflowDescriptor: EngineWorkflowDescriptor) extends ShadowWorkflowActorData

  /**
    * Mode in which the workflow should be started:
    */
  sealed trait StartMode
  case object StartNewWorkflow extends StartMode
  case object RestartExistingWorkflow extends StartMode

  def props(workflowId: WorkflowId, startMode: StartMode, wdlSource: WorkflowSourceFiles, conf: Config): Props = Props(new ShadowWorkflowActor(workflowId, startMode, wdlSource, conf))
}

/**
  * Class that orchestrates a single workflow... Shadow style.
  */
class ShadowWorkflowActor(workflowId: WorkflowId,
                          startMode: StartMode,
                          workflowSources: WorkflowSourceFiles,
                          conf: Config) extends LoggingFSM[ShadowWorkflowActorState, ShadowWorkflowActorData] with ActorLogging{

  val tag = self.path.name

  startWith(WorkflowUnstartedState, EmptyShadowWorkflowActorData)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(ShadowMaterializeWorkflowDescriptorActor.props())
      actor ! ShadowMaterializeWorkflowDescriptorCommand(workflowId, workflowSources, conf)
      goto(MaterializingWorkflowDescriptorState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(ShadowMaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), stateData) =>
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(workflowId, workflowDescriptor), name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using ShadowWorkflowActorDataWithDescriptor(workflowDescriptor)
    case Event(ShadowMaterializeWorkflowDescriptorFailureResponse(reason: Throwable), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, MaterializingWorkflowDescriptorState, Seq(reason))
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse, ShadowWorkflowActorDataWithDescriptor(workflowDescriptor)) =>
      val executionActor = context.actorOf(WorkflowExecutionActor.props(workflowId, workflowDescriptor), name = s"WorkflowExecutionActor-$workflowId")
      val commandToSend = startMode match {
        case StartNewWorkflow => StartExecutingWorkflowCommand
        case RestartExistingWorkflow => RestartExecutingWorkflowCommand
      }
      executionActor ! commandToSend
      goto(ExecutingWorkflowState)
    case Event(WorkflowInitializationFailedResponse(reason), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, InitializingWorkflowState, reason)
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse, ShadowWorkflowActorDataWithDescriptor(workflowDescriptor)) =>
      val finalizationActor = context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor), name = s"WorkflowFinalizationActor-$workflowId")
      finalizationActor ! StartEngineFinalizationCommand
      goto(FinalizingWorkflowState)
    case Event(WorkflowExecutionFailedResponse(reasons), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, ExecutingWorkflowState, reasons)
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, stateData) =>
      context.parent ! ShadowWorkflowSucceededResponse
      goto(WorkflowSucceededState)
    case Event(WorkflowFinalizationFailedResponse(reasons), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, FinalizingWorkflowState, reasons)
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  // Let these messages fall through to the whenUnhandled handler:
  when(WorkflowAbortedState) { FSM.NullFunction }
  when(WorkflowFailedState) { FSM.NullFunction }
  when(WorkflowSucceededState) { FSM.NullFunction }

  onTransition {
    case oldState -> terminalState if terminalState.terminal =>
      log.info(s"$tag transition from $oldState to $terminalState: shutting down")
      context.stop(self)
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState")
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message $unhandledMessage in state $stateName")
      stay
  }
}
