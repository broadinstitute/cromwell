package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.core.{WorkflowId, _}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine._
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationFailedResponse, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{StartInitializationCommand, WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse}
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.services.MetadataServiceActor._
import cromwell.services._

import scala.language.postfixOps

object WorkflowActor {

  /**
    * Commands to the WorkflowActor
    */
  sealed trait WorkflowActorCommand
  case object StartWorkflowCommand extends WorkflowActorCommand
  case object AbortWorkflowCommand extends WorkflowActorCommand

  /**
    * Responses from the WorkflowActor
    */
  sealed trait WorkflowActorResponse
  case class WorkflowSucceededResponse(workflowId: WorkflowId) extends WorkflowActorResponse
  case class WorkflowAbortedResponse(workflowId: WorkflowId) extends WorkflowActorResponse
  case class WorkflowFailedResponse(workflowId: WorkflowId, inState: WorkflowActorState, reasons: Seq[Throwable]) extends Exception with WorkflowActorResponse

  /**
    * States for the WorkflowActor FSM
    */
  sealed trait WorkflowActorState {
    def terminal = false
    // Each state in the FSM maps to a state in the `WorkflowState` which is used for Metadata reporting purposes
    def workflowState: WorkflowState
  }
  sealed trait WorkflowActorTerminalState extends WorkflowActorState { override val terminal = true }
  sealed trait WorkflowActorRunningState extends WorkflowActorState { override val workflowState = WorkflowRunning }

  /**
    * Waiting for a Start or Restart command.
    */
  case object WorkflowUnstartedState extends WorkflowActorState {
    override val workflowState = WorkflowSubmitted
  }

  /**
    * The WorkflowActor is created. It now needs to validate and create its WorkflowDescriptor
    */
  case object MaterializingWorkflowDescriptorState extends WorkflowActorRunningState

  /**
    * The WorkflowActor has a descriptor. It now needs to create its set of WorkflowBackendActors
    */
  case object InitializingWorkflowState extends WorkflowActorRunningState

  /**
    * The WorkflowActor has a descriptor and a set of backends. It now needs to execute the hell out of its jobs.
    */
  case object ExecutingWorkflowState extends WorkflowActorRunningState

  /**
    * The WorkflowActor has completed. So we're now finalizing whatever needs to be finalized on the backends
    */
  case object FinalizingWorkflowState extends WorkflowActorRunningState

  /**
    * The WorkflowActor is aborting. We're just waiting for everything to finish up then we'll respond back to the
    * manager.
    */
  case object WorkflowAbortingState extends WorkflowActorState {
    override val workflowState = WorkflowAborting
  }

  /**
    * We're done.
    */
  case object WorkflowSucceededState extends WorkflowActorTerminalState {
    override val workflowState = WorkflowSucceeded
  }

  /**
    * We're done. But we were aborted in the middle of the lifecycle.
    */
  case object WorkflowAbortedState extends WorkflowActorTerminalState {
    override val workflowState = WorkflowAborted
  }

  /**
    * We're done. But the workflow has failed.
    */
  case object WorkflowFailedState extends WorkflowActorTerminalState {
    override val workflowState = WorkflowFailed
  }

  case class StateCheckpoint(state: WorkflowActorState, failures: Option[List[Throwable]] = None)

  /**
    * @param currentLifecycleStateActor The current lifecycle stage, represented by an ActorRef.
    */
  case class WorkflowActorData(currentLifecycleStateActor: Option[ActorRef],
                               workflowDescriptor: Option[EngineWorkflowDescriptor],
                               lastStateReached: StateCheckpoint)
  object WorkflowActorData {
    def empty = WorkflowActorData(currentLifecycleStateActor = None, workflowDescriptor = None, StateCheckpoint(WorkflowUnstartedState))
  }

  /**
    * Mode in which the workflow should be started:
    */
  sealed trait StartMode
  case object StartNewWorkflow extends StartMode
  case object RestartExistingWorkflow extends StartMode

  def props(workflowId: WorkflowId,
            startMode: StartMode,
            wdlSource: WorkflowSourceFiles,
            conf: Config,
            serviceRegistryActor: ActorRef): Props = Props(new WorkflowActor(workflowId, startMode, wdlSource, conf, serviceRegistryActor))
}

/**
  * Class that orchestrates a single workflow.
  */
class WorkflowActor(workflowId: WorkflowId,
                    startMode: StartMode,
                    workflowSources: WorkflowSourceFiles,
                    conf: Config,
                    serviceRegistryActor: ActorRef)
  extends LoggingFSM[WorkflowActorState, WorkflowActorData] with ActorLogging {

  val tag = s"WorkflowActor [UUID(${workflowId.shortString})]"

  implicit val actorSystem = context.system
  implicit val ec = context.dispatcher

  startWith(WorkflowUnstartedState, WorkflowActorData.empty)

  pushCurrentStateToMetadataService(WorkflowUnstartedState.workflowState)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(MaterializeWorkflowDescriptorActor.props(), s"MaterializeWorkflowDescriptorActor-$workflowId")
      actor ! MaterializeWorkflowDescriptorCommand(workflowId, workflowSources, conf)
      goto(MaterializingWorkflowDescriptorState) using stateData.copy(currentLifecycleStateActor = Option(actor))
    case Event(AbortWorkflowCommand, stateData) => goto(WorkflowAbortedState)
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), data) =>
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(workflowId, workflowDescriptor), name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using data.copy(currentLifecycleStateActor = Option(initializerActor), workflowDescriptor = Option(workflowDescriptor))
    case Event(MaterializeWorkflowDescriptorFailureResponse(reason: Throwable), _) =>
      context.parent ! WorkflowFailedResponse(workflowId, MaterializingWorkflowDescriptorState, Seq(reason))
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) =>
      // No lifecycle sub-actors exist yet, so no indirection via WorkflowAbortingState is necessary:
      sender ! WorkflowAbortedResponse(workflowId)
      goto(WorkflowAbortedState)
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse, data @ WorkflowActorData(_, Some(workflowDescriptor), _)) =>
      // Add Workflow name and inputs to the metadata
      pushWfNameAndInputsToMetadataService(workflowDescriptor)

      val executionActor = context.actorOf(WorkflowExecutionActor.props(workflowId, workflowDescriptor, serviceRegistryActor), name = s"WorkflowExecutionActor-$workflowId")
      val commandToSend = startMode match {
        case StartNewWorkflow => StartExecutingWorkflowCommand
        case RestartExistingWorkflow => RestartExecutingWorkflowCommand
      }
      executionActor ! commandToSend
      goto(ExecutingWorkflowState) using data.copy(currentLifecycleStateActor = Option(executionActor))
    case Event(WorkflowInitializationFailedResponse(reason), data @ WorkflowActorData(_, Some(workflowDescriptor), _)) =>
      finalizeWorkflow(data, workflowDescriptor, failures = Option(reason.toList))
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse, data @ WorkflowActorData(_, Some(workflowDescriptor), _)) =>
      finalizeWorkflow(data, workflowDescriptor, failures = None)
    case Event(WorkflowExecutionFailedResponse(failures), data @ WorkflowActorData(_, Some(workflowDescriptor), _)) =>
      finalizeWorkflow(data, workflowDescriptor, failures = Option(failures.toList))
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, data) => finalizationSucceeded(data)
    case Event(WorkflowFinalizationFailedResponse(finalizationFailures), data) =>
      val failures = data.lastStateReached.failures.getOrElse(List.empty) ++ finalizationFailures
      context.parent ! WorkflowFailedResponse(workflowId, FinalizingWorkflowState, failures)

      goto(WorkflowFailedState)
  }

  when(WorkflowAbortingState) {
    case Event(x: EngineLifecycleStateCompleteResponse, data @ WorkflowActorData(_, Some(workflowDescriptor), _)) =>
      finalizeWorkflow(data, workflowDescriptor, failures = None)
    case _ => stay()
  }

  // Let these messages fall through to the whenUnhandled handler:
  when(WorkflowAbortedState) { FSM.NullFunction }
  when(WorkflowFailedState) { FSM.NullFunction }
  when(WorkflowSucceededState) { FSM.NullFunction }

  whenUnhandled {
    case Event(MetadataPutFailed(action, error), _) =>
      // Do something useful here??
      log.warning(s"$tag Put failed for Metadata action $action : ${error.getMessage}")
      stay
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    case Event(AbortWorkflowCommand, WorkflowActorData(Some(actor), _, _)) =>
      actor ! EngineLifecycleActorAbortCommand
      goto(WorkflowAbortingState)
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message $unhandledMessage in state $stateName")
      stay
  }

  onTransition {
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState")
      // This updates the workflow status
      // Only publish "External" state to metadata service
      // workflowState maps a state to an "external" state (e.g all states extending WorkflowActorRunningState map to WorkflowRunning)
      if (fromState.workflowState != toState.workflowState) {
        pushCurrentStateToMetadataService(toState.workflowState)
      }
  }

  onTransition {
    case oldState -> terminalState if terminalState.terminal =>
      log.info(s"$tag transition from $oldState to $terminalState: shutting down")
      // Add the end time of the workflow in the MetadataService
      val now = OffsetDateTime.now
      val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime),
        MetadataValue(now), now)
      serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
      terminalState match {
        case WorkflowSucceededState =>
          context.parent ! WorkflowSucceededResponse(workflowId)
        case WorkflowFailedState => // The failure will already have been sent. So just kick back and relax here.
        case WorkflowAbortedState =>
          context.parent ! WorkflowAbortedResponse(workflowId)
        case unknownState => log.warning(s"$tag $unknownState is an unhandled terminal state!")
      }
  }

  private def finalizationSucceeded(data: WorkflowActorData) = {
    val finalState = data.lastStateReached match {
      case StateCheckpoint(WorkflowAbortingState, None) =>
        context.parent ! WorkflowAbortedResponse(workflowId)
        WorkflowAbortedState
      case StateCheckpoint(state, Some(failures)) =>
        context.parent ! WorkflowFailedResponse(workflowId, state, failures)
        WorkflowFailedState
      case StateCheckpoint(state, None) =>
        context.parent ! WorkflowSucceededResponse(workflowId)
        WorkflowSucceededState
    }

    goto(finalState) using data.copy(currentLifecycleStateActor = None)
  }

  /**
    * Run finalization actor and transition to FinalizingWorkflowState.
    */
  private def finalizeWorkflow(data: WorkflowActorData, workflowDescriptor: EngineWorkflowDescriptor, failures: Option[List[Throwable]]) = {
    val finalizationActor = context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor), name = s"WorkflowFinalizationActor-$workflowId")
    finalizationActor ! StartFinalizationCommand
    goto(FinalizingWorkflowState) using data.copy(lastStateReached = StateCheckpoint(stateName, failures))
  }

  private def pushWfNameAndInputsToMetadataService(workflowDescriptor: EngineWorkflowDescriptor): Unit = {
    import MetadataServiceActorImplicits.EnhancedServiceRegistryActorForMetadata
    // Inputs
    workflowDescriptor.backendDescriptor.inputs.foreach { case (inputName, wdlValue) =>
      serviceRegistryActor.pushWdlValueMetadata(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Inputs}:$inputName"), wdlValue)
    }
    // Workflow name:
    val nameEvent = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Name), MetadataValue(workflowDescriptor.name), OffsetDateTime.now)
    serviceRegistryActor ! PutMetadataAction(nameEvent)
  }

  // Update the current State of the Workflow (corresponding to the FSM state) in the Metadata service
  private def pushCurrentStateToMetadataService(workflowState: WorkflowState): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status),
      MetadataValue(workflowState), OffsetDateTime.now)
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }
}
