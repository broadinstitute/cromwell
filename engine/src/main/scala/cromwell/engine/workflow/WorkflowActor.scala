package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.backend.AllBackendInitializationData
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.FinalWorkflowLogDir
import cromwell.core._
import cromwell.core.logging.{WorkflowLogger, WorkflowLogging}
import cromwell.core.path.PathFactory
import cromwell.engine._
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationFailedResponse, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{StartInitializationCommand, WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse}
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.webservice.EngineStatsActor

import scala.util.Random

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
                               initializationData: AllBackendInitializationData,
                               lastStateReached: StateCheckpoint)
  object WorkflowActorData {
    def empty = WorkflowActorData(
      currentLifecycleStateActor = None,
      workflowDescriptor = None,
      initializationData = AllBackendInitializationData.empty,
      lastStateReached = StateCheckpoint(WorkflowUnstartedState))
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
            serviceRegistryActor: ActorRef,
            workflowLogCopyRouter: ActorRef,
            jobStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection): Props = {
    Props(new WorkflowActor(workflowId, startMode, wdlSource, conf, serviceRegistryActor, workflowLogCopyRouter,
      jobStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection)).withDispatcher(EngineDispatcher)
  }
}

/**
  * Class that orchestrates a single workflow.
  */
class WorkflowActor(val workflowId: WorkflowId,
                    startMode: StartMode,
                    workflowSources: WorkflowSourceFiles,
                    conf: Config,
                    serviceRegistryActor: ActorRef,
                    workflowLogCopyRouter: ActorRef,
                    jobStoreActor: ActorRef,
                    callCacheReadActor: ActorRef,
                    jobTokenDispenserActor: ActorRef,
                    backendSingletonCollection: BackendSingletonCollection)
  extends LoggingFSM[WorkflowActorState, WorkflowActorData] with WorkflowLogging {

  implicit val ec = context.dispatcher

  startWith(WorkflowUnstartedState, WorkflowActorData.empty)

  pushCurrentStateToMetadataService(WorkflowUnstartedState.workflowState)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryActor, workflowId),
        "MaterializeWorkflowDescriptorActor")
      val startEvent = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.StartTime), MetadataValue(OffsetDateTime.now.toString))
      serviceRegistryActor ! PutMetadataAction(startEvent)

      actor ! MaterializeWorkflowDescriptorCommand(workflowSources, conf)
      goto(MaterializingWorkflowDescriptorState) using stateData.copy(currentLifecycleStateActor = Option(actor))
    case Event(AbortWorkflowCommand, stateData) => goto(WorkflowAbortedState)
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), data) =>
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(workflowId, workflowDescriptor, serviceRegistryActor),
        name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using data.copy(currentLifecycleStateActor = Option(initializerActor), workflowDescriptor = Option(workflowDescriptor))
    case Event(MaterializeWorkflowDescriptorFailureResponse(reason: Throwable), data) =>
      goto(WorkflowFailedState) using data.copy(lastStateReached = StateCheckpoint(MaterializingWorkflowDescriptorState, Option(List(reason))))
    case Event(AbortWorkflowCommand, stateData) =>
      // No lifecycle sub-actors exist yet, so no indirection via WorkflowAbortingState is necessary:
      goto(WorkflowAbortedState)
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse(initializationData), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      val restarting = startMode match {
        case StartNewWorkflow => false
        case RestartExistingWorkflow => true
      }

      val executionActor = context.actorOf(WorkflowExecutionActor.props(workflowId,
        workflowDescriptor,
        serviceRegistryActor,
        jobStoreActor,
        callCacheReadActor,
        jobTokenDispenserActor,
        backendSingletonCollection,
        initializationData,
        restarting = restarting), name = s"WorkflowExecutionActor-$workflowId")

      executionActor ! ExecuteWorkflowCommand

      goto(ExecutingWorkflowState) using data.copy(currentLifecycleStateActor = Option(executionActor), initializationData = initializationData)
    case Event(WorkflowInitializationFailedResponse(reason), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, ExecutionStore.empty, OutputStore.empty, Option(reason.toList))
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse(executionStore, outputStore),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, executionStore, outputStore, None)
    case Event(WorkflowExecutionFailedResponse(executionStore, outputStore, failures),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, executionStore, outputStore, Option(failures.toList))
    case Event(msg @ EngineStatsActor.JobCountQuery, data) =>
      data.currentLifecycleStateActor match {
        case Some(a) => a forward msg
        case None => sender ! EngineStatsActor.NoJobs // This should be impossible, but if somehow here it's technically correct
      }

      stay()
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, data) => finalizationSucceeded(data)
    case Event(WorkflowFinalizationFailedResponse(finalizationFailures), data) =>
      val failures = data.lastStateReached.failures.getOrElse(List.empty) ++ finalizationFailures
      goto(WorkflowFailedState) using data.copy(lastStateReached = StateCheckpoint(FinalizingWorkflowState, Option(failures)))
  }

  when(WorkflowAbortingState) {
    case Event(x: EngineLifecycleStateCompleteResponse, data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, ExecutionStore.empty, OutputStore.empty, failures = None)
    case _ => stay()
  }

  // Let these messages fall through to the whenUnhandled handler:
  when(WorkflowAbortedState) { FSM.NullFunction }
  when(WorkflowFailedState) { FSM.NullFunction }
  when(WorkflowSucceededState) { FSM.NullFunction }

  whenUnhandled {
    case Event(MetadataPutFailed(action, error), _) =>
      // Do something useful here??
      workflowLogger.warn(s"Put failed for Metadata action $action : ${error.getMessage}")
      stay
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    case Event(AbortWorkflowCommand, WorkflowActorData(Some(actor), _, _, _)) =>
      actor ! EngineLifecycleActorAbortCommand
      goto(WorkflowAbortingState)
    case Event(EngineStatsActor.JobCountQuery, _) =>
      sender ! EngineStatsActor.NoJobs
      stay()
    case unhandledMessage =>
      workflowLogger.warn(s"received an unhandled message $unhandledMessage in state $stateName")
      stay
  }

  onTransition {
    case fromState -> toState =>
      workflowLogger.debug(s"transitioning from {} to {}", arg1 = fromState, arg2 = toState)
      // This updates the workflow status
      // Only publish "External" state to metadata service
      // workflowState maps a state to an "external" state (e.g all states extending WorkflowActorRunningState map to WorkflowRunning)
      if (fromState.workflowState != toState.workflowState) {
        pushCurrentStateToMetadataService(toState.workflowState)
      }
  }

  onTransition {
    case (oldState, terminalState: WorkflowActorTerminalState) =>
      workflowLogger.debug(s"transition from {} to {}. Stopping self.", arg1 = oldState, arg2 = terminalState)
      // Add the end time of the workflow in the MetadataService
      val now = OffsetDateTime.now
      val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime), MetadataValue(now))
      serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
      terminalState match {
        case WorkflowFailedState =>
          val failures = nextStateData.lastStateReached.failures.getOrElse(List.empty)
          val failureEvents = failures flatMap { r => throwableToMetadataEvents(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Failures}[${Random.nextInt(Int.MaxValue)}]"), r) }
          serviceRegistryActor ! PutMetadataAction(failureEvents)
          context.parent ! WorkflowFailedResponse(workflowId, nextStateData.lastStateReached.state, failures)
        case _ => // The WMA is watching state transitions and needs no further info
      }

      // Copy/Delete workflow logs
      if (WorkflowLogger.isEnabled) {
        stateData.workflowDescriptor foreach { wd =>
          wd.getWorkflowOption(FinalWorkflowLogDir) match {
            case Some(destinationDir) =>
              workflowLogCopyRouter ! CopyWorkflowLogsActor.Copy(wd.id, PathFactory.buildPath(destinationDir, wd.pathBuilders))
            case None if WorkflowLogger.isTemporary => workflowLogger.deleteLogFile()
            case _ =>
          }
        }
      }

      context stop self
  }

  private def finalizationSucceeded(data: WorkflowActorData) = {
    val finalState = data.lastStateReached match {
      case StateCheckpoint(WorkflowAbortingState, None) => WorkflowAbortedState
      case StateCheckpoint(state, Some(failures)) => WorkflowFailedState
      case StateCheckpoint(state, None) => WorkflowSucceededState
    }

    goto(finalState) using data.copy(currentLifecycleStateActor = None)
  }

  private[workflow] def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, executionStore: ExecutionStore, outputStore: OutputStore) = {
    context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor, executionStore, outputStore, stateData.initializationData), name = s"WorkflowFinalizationActor")
  }
  /**
    * Run finalization actor and transition to FinalizingWorkflowState.
    */
  private def finalizeWorkflow(data: WorkflowActorData, workflowDescriptor: EngineWorkflowDescriptor,
                               executionStore: ExecutionStore, outputStore: OutputStore,
                               failures: Option[List[Throwable]]) = {
    val finalizationActor = makeFinalizationActor(workflowDescriptor, executionStore, outputStore)
    finalizationActor ! StartFinalizationCommand
    goto(FinalizingWorkflowState) using data.copy(lastStateReached = StateCheckpoint(stateName, failures))
  }

  // Update the current State of the Workflow (corresponding to the FSM state) in the Metadata service
  private def pushCurrentStateToMetadataService(workflowState: WorkflowState): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status),
      MetadataValue(workflowState))
    serviceRegistryActor ! PutMetadataAction(metadataEventMsg)
  }
}
