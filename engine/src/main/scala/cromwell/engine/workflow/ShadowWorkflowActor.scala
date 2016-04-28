package cromwell.engine.workflow

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.{EngineWorkflowDescriptor, WorkflowSourceFiles}

import cromwell.engine.workflow.lifecycle.ShadowMaterializeWorkflowDescriptorActor.{ShadowMaterializeWorkflowDescriptorFailureResponse, ShadowMaterializeWorkflowDescriptorSuccessResponse, ShadowMaterializeWorkflowDescriptorCommand}
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.{WorkflowExecutionSucceededResponse, WorkflowExecutionFailedResponse, RestartExecutingWorkflowCommand, StartExecutingWorkflowCommand}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{WorkflowFinalizationSucceededResponse, WorkflowFinalizationFailedResponse, StartFinalizationCommand}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse, StartInitializationCommand}
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.ShadowWorkflowActor._
import scala.language.postfixOps

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
  case class ShadowWorkflowSucceededResponse(workflowId: WorkflowId) extends ShadowWorkflowActorResponse
  case class ShadowWorkflowAbortedResponse(workflowId: WorkflowId) extends ShadowWorkflowActorResponse
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
    * The WorkflowActor is aborting. We're just waiting for everything to finish up then we'll respond back to the
    * manager.
    */
  case object AbortingWorkflowState extends ShadowWorkflowActorState

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

  /**
    * @param currentLifecycleStateActor The current lifecycle stage, represented by an ActorRef.
    */
  case class ShadowWorkflowActorData(currentLifecycleStateActor: Option[ActorRef],
                                     workflowDescriptor: Option[EngineWorkflowDescriptor])
  object ShadowWorkflowActorData {
    def empty = ShadowWorkflowActorData(None, None)
    def apply(currentLifecycleStateActor: ActorRef, workflowDescriptor: EngineWorkflowDescriptor): ShadowWorkflowActorData = ShadowWorkflowActorData(Option(currentLifecycleStateActor), Option(workflowDescriptor))
  }

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

  startWith(WorkflowUnstartedState, ShadowWorkflowActorData.empty)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(ShadowMaterializeWorkflowDescriptorActor.props())
      actor ! ShadowMaterializeWorkflowDescriptorCommand(workflowId, workflowSources, conf)
      goto(MaterializingWorkflowDescriptorState) using stateData.copy(currentLifecycleStateActor = Option(actor))
    case Event(AbortWorkflowCommand, stateData) =>
      // No lifecycle sub-actors exist yet, so no indirection via WorkflowAbortingState is necessary:
      sender ! ShadowWorkflowAbortedResponse(workflowId)
      goto(WorkflowAbortedState)
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(ShadowMaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), stateData) =>
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(workflowId, workflowDescriptor), name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using ShadowWorkflowActorData(initializerActor, workflowDescriptor)
    case Event(ShadowMaterializeWorkflowDescriptorFailureResponse(reason: Throwable), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, MaterializingWorkflowDescriptorState, Seq(reason))
      goto(WorkflowFailedState)
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse, ShadowWorkflowActorData(_, Some(workflowDescriptor))) =>
      val executionActor = context.actorOf(WorkflowExecutionActor.props(workflowId, workflowDescriptor), name = s"WorkflowExecutionActor-$workflowId")
      val commandToSend = startMode match {
        case StartNewWorkflow => StartExecutingWorkflowCommand
        case RestartExistingWorkflow => RestartExecutingWorkflowCommand
      }
      executionActor ! commandToSend
      goto(ExecutingWorkflowState) using ShadowWorkflowActorData(executionActor, workflowDescriptor)
    case Event(WorkflowInitializationFailedResponse(reason), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, InitializingWorkflowState, reason)
      goto(WorkflowFailedState)
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse, ShadowWorkflowActorData(_, Some(workflowDescriptor))) =>
      val finalizationActor = context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor), name = s"WorkflowFinalizationActor-$workflowId")
      finalizationActor ! StartFinalizationCommand
      goto(FinalizingWorkflowState) using ShadowWorkflowActorData(finalizationActor, workflowDescriptor)
    case Event(WorkflowExecutionFailedResponse(reasons), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, ExecutingWorkflowState, reasons)
      goto(WorkflowFailedState)
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, stateData) =>
      context.parent ! ShadowWorkflowSucceededResponse(workflowId)
      goto(WorkflowSucceededState)
    case Event(WorkflowFinalizationFailedResponse(reasons), _) =>
      context.parent ! ShadowWorkflowFailedResponse(workflowId, FinalizingWorkflowState, reasons)
      goto(WorkflowFailedState)
  }

  when(AbortingWorkflowState) {
    case Event(x: EngineLifecycleStateCompleteResponse, _) =>
      context.parent ! ShadowWorkflowAbortedResponse(workflowId)
      goto(WorkflowAbortedState)
    case _ => stay()
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
    case Event(AbortWorkflowCommand, ShadowWorkflowActorData(Some(actor), _)) =>
      actor ! EngineLifecycleActorAbortCommand
      goto(AbortingWorkflowState)
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message $unhandledMessage in state $stateName")
      stay
  }
}
