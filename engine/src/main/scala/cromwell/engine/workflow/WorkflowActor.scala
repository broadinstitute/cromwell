package cromwell.engine.workflow

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.core.{KnowsWhatTimeItIs, WorkflowId}
import cromwell.engine._
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationFailedResponse, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{StartInitializationCommand, WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse}
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{RestartExecutingWorkflowCommand, StartExecutingWorkflowCommand, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.services.MetadataServiceActor._
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue}
import org.joda.time.DateTime

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
  case object AbortingWorkflowState extends WorkflowActorState {
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

  /**
    * @param currentLifecycleStateActor The current lifecycle stage, represented by an ActorRef.
    */
  case class WorkflowActorData(currentLifecycleStateActor: Option[ActorRef],
                               workflowDescriptor: Option[EngineWorkflowDescriptor])
  object WorkflowActorData {
    def empty = WorkflowActorData(None, None)
    def apply(currentLifecycleStateActor: ActorRef, workflowDescriptor: EngineWorkflowDescriptor): WorkflowActorData = WorkflowActorData(Option(currentLifecycleStateActor), Option(workflowDescriptor))
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
  extends LoggingFSM[WorkflowActorState, WorkflowActorData] with ActorLogging with KnowsWhatTimeItIs {

  val tag = self.path.name

  implicit val actorSystem = context.system
  implicit val ec = context.dispatcher

  startWith(WorkflowUnstartedState, WorkflowActorData.empty)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(MaterializeWorkflowDescriptorActor.props())
      actor ! MaterializeWorkflowDescriptorCommand(workflowId, workflowSources, conf)
      goto(MaterializingWorkflowDescriptorState) using stateData.copy(currentLifecycleStateActor = Option(actor))
    case Event(AbortWorkflowCommand, stateData) =>
      // No lifecycle sub-actors exist yet, so no indirection via WorkflowAbortingState is necessary:
      sender ! WorkflowAbortedResponse(workflowId)
      goto(WorkflowAbortedState)
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), stateData) =>
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(workflowId, workflowDescriptor), name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using WorkflowActorData(initializerActor, workflowDescriptor)
    case Event(MaterializeWorkflowDescriptorFailureResponse(reason: Throwable), _) =>
      context.parent ! WorkflowFailedResponse(workflowId, MaterializingWorkflowDescriptorState, Seq(reason))
      goto(WorkflowFailedState)
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse, WorkflowActorData(_, Some(workflowDescriptor))) =>
      // Add Workflow name and inputs to the metadata
      pushWfNameAndInputsToMetadataService(workflowDescriptor)

      val executionActor = context.actorOf(WorkflowExecutionActor.props(workflowId, workflowDescriptor, serviceRegistryActor), name = s"WorkflowExecutionActor-$workflowId")
      val commandToSend = startMode match {
        case StartNewWorkflow => StartExecutingWorkflowCommand
        case RestartExistingWorkflow => RestartExecutingWorkflowCommand
      }
      executionActor ! commandToSend
      goto(ExecutingWorkflowState) using WorkflowActorData(executionActor, workflowDescriptor)
    case Event(WorkflowInitializationFailedResponse(reason), _) =>
      context.parent ! WorkflowFailedResponse(workflowId, InitializingWorkflowState, reason)
      goto(WorkflowFailedState)
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse, WorkflowActorData(_, Some(workflowDescriptor))) =>
      val finalizationActor = context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor), name = s"WorkflowFinalizationActor-$workflowId")
      finalizationActor ! StartFinalizationCommand
      goto(FinalizingWorkflowState) using WorkflowActorData(finalizationActor, workflowDescriptor)
    case Event(WorkflowExecutionFailedResponse(reasons), _) =>
      context.parent ! WorkflowFailedResponse(workflowId, ExecutingWorkflowState, reasons)
      goto(WorkflowFailedState)
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, stateData) =>
      context.parent ! WorkflowSucceededResponse(workflowId)
      goto(WorkflowSucceededState)
    case Event(WorkflowFinalizationFailedResponse(reasons), _) =>
      context.parent ! WorkflowFailedResponse(workflowId, FinalizingWorkflowState, reasons)
      goto(WorkflowFailedState)
  }

  when(AbortingWorkflowState) {
    case Event(x: EngineLifecycleStateCompleteResponse, _) =>
      context.parent ! WorkflowAbortedResponse(workflowId)
      goto(WorkflowAbortedState)
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
    case Event(AbortWorkflowCommand, WorkflowActorData(Some(actor), _)) =>
      actor ! EngineLifecycleActorAbortCommand
      goto(AbortingWorkflowState)
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message $unhandledMessage in state $stateName")
      stay
  }

  onTransition {
    // Only publish "External" state to metadata service
    // workflowState maps a state to an "external" state (e.g all states extending WorkflowActorRunningState map to WorkflowRunning)
    case fromState -> toState if fromState.workflowState != toState.workflowState =>
      log.info(s"$tag transitioning from $fromState to $toState")
      // This updates the workflow status
      pushCurrentStateToMetadataService(toState.workflowState)
  }

  onTransition {
    case oldState -> terminalState if terminalState.terminal =>
      log.info(s"$tag transition from $oldState to $terminalState: shutting down")
      // Add the end time of the workflow in the MetadataService
      import KnowsWhatTimeItIs._
      val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime), MetadataValue(currentTime.asJodaString), currentTime)
      serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
      terminalState match {
        case WorkflowSucceededState =>
          context.parent ! WorkflowSucceededResponse(workflowId)
        case WorkflowFailedState =>
          context.parent ! WorkflowFailedResponse(workflowId, oldState, Seq.empty[Throwable])
        case WorkflowAbortedState =>
          context.parent ! WorkflowAbortedResponse(workflowId)
        case unknownState => log.warning(s"$tag $unknownState is an unhandled terminal state!")
      }
  }

  private def pushWfNameAndInputsToMetadataService(workflowDescriptor: EngineWorkflowDescriptor): Unit = {
    val inputMetadataEvents = workflowDescriptor.backendDescriptor.inputs.map { case (k, v) =>
      MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Inputs}:$k"), MetadataValue(v.toWdlString), currentTime)
    }
    val metadataEventMsgs = List(MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Name), MetadataValue(workflowDescriptor.name), currentTime)) ++ inputMetadataEvents
    metadataEventMsgs foreach ( serviceRegistryActor ! PutMetadataAction(_) )
  }

  // Update the current State of the Workflow (corresponding to the FSM state) in the Metadata service
  private def pushCurrentStateToMetadataService(workflowState: WorkflowState): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue(workflowState.toString), currentTime)
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }

}
