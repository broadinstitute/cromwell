package cromwell.engine.workflow

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.FinalWorkflowLogDir
import cromwell.core._
import cromwell.core.logging.{WorkflowLogger, WorkflowLogging}
import cromwell.core.path.{PathBuilder, PathFactory}
import cromwell.engine._
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationFailedResponse, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{StartInitializationCommand, WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse}
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.lifecycle.execution.{WorkflowExecutionActor, WorkflowMetadataHelper}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.services.metadata.MetadataService._
import cromwell.subworkflowstore.SubWorkflowStoreActor.WorkflowComplete
import cromwell.webservice.EngineStatsActor
import wdl4s.{LocallyQualifiedName => _}

import scala.util.Failure

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
            wdlSource: WorkflowSourceFilesCollection,
            conf: Config,
            serviceRegistryActor: ActorRef,
            workflowLogCopyRouter: ActorRef,
            jobStoreActor: ActorRef,
            subWorkflowStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection,
            serverMode: Boolean): Props = {
    Props(new WorkflowActor(workflowId, startMode, wdlSource, conf, serviceRegistryActor, workflowLogCopyRouter,
      jobStoreActor, subWorkflowStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection, serverMode)).withDispatcher(EngineDispatcher)
  }
}

/**
  * Class that orchestrates a single workflow.
  */
class WorkflowActor(val workflowId: WorkflowId,
                    startMode: StartMode,
                    workflowSources: WorkflowSourceFilesCollection,
                    conf: Config,
                    override val serviceRegistryActor: ActorRef,
                    workflowLogCopyRouter: ActorRef,
                    jobStoreActor: ActorRef,
                    subWorkflowStoreActor: ActorRef,
                    callCacheReadActor: ActorRef,
                    jobTokenDispenserActor: ActorRef,
                    backendSingletonCollection: BackendSingletonCollection,
                    serverMode: Boolean)
  extends LoggingFSM[WorkflowActorState, WorkflowActorData] with WorkflowLogging with WorkflowMetadataHelper {

  implicit val ec = context.dispatcher
  override val workflowIdForLogging = workflowId

  startWith(WorkflowUnstartedState, WorkflowActorData.empty)

  pushCurrentStateToMetadataService(workflowId, WorkflowUnstartedState.workflowState)
  
  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryActor, workflowId, importLocalFilesystem = !serverMode),
        "MaterializeWorkflowDescriptorActor")
      pushWorkflowStart(workflowId)
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

      val executionActor = context.actorOf(WorkflowExecutionActor.props(
        workflowDescriptor,
        serviceRegistryActor,
        jobStoreActor,
        subWorkflowStoreActor,
        callCacheReadActor,
        jobTokenDispenserActor,
        backendSingletonCollection,
        initializationData,
        restarting = restarting), name = s"WorkflowExecutionActor-$workflowId")

      executionActor ! ExecuteWorkflowCommand

      goto(ExecutingWorkflowState) using data.copy(currentLifecycleStateActor = Option(executionActor), initializationData = initializationData)
    case Event(WorkflowInitializationFailedResponse(reason), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, Map.empty, Map.empty, Option(reason.toList))
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionSucceededResponse(jobKeys, outputs),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, jobKeys, outputs, None)
    case Event(WorkflowExecutionFailedResponse(jobKeys, failures),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, jobKeys, Map.empty, Option(List(failures)))
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
      finalizeWorkflow(data, workflowDescriptor, Map.empty, Map.empty, failures = None)
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
        pushCurrentStateToMetadataService(workflowId, toState.workflowState)
      }
  }

  onTransition {
    case (oldState, terminalState: WorkflowActorTerminalState) =>
      workflowLogger.debug(s"transition from {} to {}. Stopping self.", arg1 = oldState, arg2 = terminalState)
      pushWorkflowEnd(workflowId)
      subWorkflowStoreActor ! WorkflowComplete(workflowId)
      terminalState match {
        case WorkflowFailedState =>
          val failures = nextStateData.lastStateReached.failures.getOrElse(List.empty)
          pushWorkflowFailures(workflowId, failures)
          context.parent ! WorkflowFailedResponse(workflowId, nextStateData.lastStateReached.state, failures)
        case _ => // The WMA is watching state transitions and needs no further info
      }

      // Copy/Delete workflow logs
      if (WorkflowLogger.isEnabled) {
        /*
          * The submitted workflow options have been previously validated by the CromwellApiHandler. These are
          * being recreated so that in case MaterializeWorkflowDescriptor fails, the workflow logs can still
          * be copied by accessing the workflow options outside of the EngineWorkflowDescriptor.
          */
        def bruteForceWorkflowOptions: WorkflowOptions = WorkflowOptions.fromJsonString(workflowSources.workflowOptionsJson).getOrElse(WorkflowOptions.fromJsonString("{}").get)
        def bruteForcePathBuilders: List[PathBuilder] = EngineFilesystems(context.system).pathBuildersForWorkflow(bruteForceWorkflowOptions)

        val (workflowOptions, pathBuilders) = stateData.workflowDescriptor match {
          case Some(wd) => (wd.backendDescriptor.workflowOptions, wd.pathBuilders)
          case None => (bruteForceWorkflowOptions, bruteForcePathBuilders)
        }

        workflowOptions.get(FinalWorkflowLogDir).toOption match {
            case Some(destinationDir) =>
              workflowLogCopyRouter ! CopyWorkflowLogsActor.Copy(workflowId, PathFactory.buildPath(destinationDir, pathBuilders))
            case None if WorkflowLogger.isTemporary => workflowLogger.deleteLogFile() match {
              case Failure(f) => log.error(f, "Failed to delete workflow log")
              case _ =>
            }
            case _ =>
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

  private[workflow] def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs) = {
    context.actorOf(WorkflowFinalizationActor.props(workflowId, workflowDescriptor, jobExecutionMap, workflowOutputs, stateData.initializationData), name = s"WorkflowFinalizationActor")
  }
  /**
    * Run finalization actor and transition to FinalizingWorkflowState.
    */
  private def finalizeWorkflow(data: WorkflowActorData, workflowDescriptor: EngineWorkflowDescriptor,
                               jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                               failures: Option[List[Throwable]]) = {
    val finalizationActor = makeFinalizationActor(workflowDescriptor, jobExecutionMap, workflowOutputs)
    finalizationActor ! StartFinalizationCommand
    goto(FinalizingWorkflowState) using data.copy(lastStateReached = StateCheckpoint (stateName, failures))
  }

}
