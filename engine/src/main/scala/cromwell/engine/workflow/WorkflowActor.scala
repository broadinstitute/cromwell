package cromwell.engine.workflow

import java.util.concurrent.atomic.AtomicInteger

import akka.actor.SupervisorStrategy.Stop
import akka.actor._
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.standard.callcaching.BlacklistCache
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowOptions.FinalWorkflowLogDir
import cromwell.core.WorkflowProcessingEvents.DescriptionEventValue.Finished
import cromwell.core._
import cromwell.core.logging.{WorkflowLogger, WorkflowLogging}
import cromwell.core.path.{Path, PathBuilder, PathBuilderFactory, PathFactory}
import cromwell.engine._
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.instrumentation.WorkflowInstrumentation
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle._
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor.{DeleteWorkflowFilesFailedResponse, DeleteWorkflowFilesSucceededResponse, StartWorkflowFilesDeletion}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.finalization.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationFailedResponse, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.finalization.{CopyWorkflowLogsActor, CopyWorkflowOutputsActor, WorkflowFinalizationActor}
import cromwell.engine.workflow.lifecycle.initialization.WorkflowInitializationActor
import cromwell.engine.workflow.lifecycle.initialization.WorkflowInitializationActor.{StartInitializationCommand, WorkflowInitializationFailedResponse, WorkflowInitializationResponse, WorkflowInitializationSucceededResponse}
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.WorkflowStoreWriteHeartbeatCommand
import cromwell.engine.workflow.workflowstore.{RestartableAborting, StartableState, WorkflowHeartbeatConfig, WorkflowToStart}
import cromwell.subworkflowstore.SubWorkflowStoreActor.WorkflowComplete
import cromwell.webservice.EngineStatsActor
import org.apache.commons.lang3.exception.ExceptionUtils
import wom.values.WomValue

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.Failure

object WorkflowActor {

  /**
    * Commands to the WorkflowActor
    */
  sealed trait WorkflowActorCommand
  case object StartWorkflowCommand extends WorkflowActorCommand
  case object AbortWorkflowCommand extends WorkflowActorCommand
  final case class AbortWorkflowWithExceptionCommand(exception: Throwable) extends WorkflowActorCommand
  case object SendWorkflowHeartbeatCommand extends WorkflowActorCommand

  case class WorkflowFailedResponse(workflowId: WorkflowId, inState: WorkflowActorState, reasons: Seq[Throwable])

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
    * The workflow and finalization has succeeded and wants to delete intermediate files. So we are now in deleting
    * those files state.
    */
  case object DeletingFilesState extends WorkflowActorRunningState

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
                               lastStateReached: StateCheckpoint,
                               effectiveStartableState: StartableState,
                               workflowFinalOutputs: Set[WomValue] = Set.empty,
                               workflowAllOutputs: Set[WomValue] = Set.empty,
                               rootAndSubworkflowIds: Set[WorkflowId] = Set.empty)
  object WorkflowActorData {
    def apply(startableState: StartableState): WorkflowActorData = WorkflowActorData(
      currentLifecycleStateActor = None,
      workflowDescriptor = None,
      initializationData = AllBackendInitializationData.empty,
      lastStateReached = StateCheckpoint(WorkflowUnstartedState),
      effectiveStartableState = startableState)
  }

  /**
    * Mode in which the workflow should be started:
    */
  sealed trait StartMode
  case object StartNewWorkflow extends StartMode
  case object RestartExistingWorkflow extends StartMode

  def props(workflowToStart: WorkflowToStart,
            conf: Config,
            callCachingEnabled: Boolean,
            invalidateBadCacheResults: Boolean,
            ioActor: ActorRef,
            serviceRegistryActor: ActorRef,
            workflowLogCopyRouter: ActorRef,
            jobStoreActor: ActorRef,
            subWorkflowStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            callCacheWriteActor: ActorRef,
            dockerHashActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            workflowStoreActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection,
            serverMode: Boolean,
            workflowHeartbeatConfig: WorkflowHeartbeatConfig,
            totalJobsByRootWf: AtomicInteger,
            fileHashCacheActorProps: Option[Props],
            blacklistCache: Option[BlacklistCache]): Props = {
    Props(
      new WorkflowActor(
        workflowToStart = workflowToStart,
        conf = conf,
        callCachingEnabled = callCachingEnabled,
        invalidateBadCacheResults = invalidateBadCacheResults,
        ioActor = ioActor,
        serviceRegistryActor = serviceRegistryActor,
        workflowLogCopyRouter = workflowLogCopyRouter,
        jobStoreActor = jobStoreActor,
        subWorkflowStoreActor = subWorkflowStoreActor,
        callCacheReadActor = callCacheReadActor,
        callCacheWriteActor = callCacheWriteActor,
        dockerHashActor = dockerHashActor,
        jobTokenDispenserActor = jobTokenDispenserActor,
        workflowStoreActor = workflowStoreActor,
        backendSingletonCollection = backendSingletonCollection,
        serverMode = serverMode,
        workflowHeartbeatConfig = workflowHeartbeatConfig,
        totalJobsByRootWf = totalJobsByRootWf,
        fileHashCacheActorProps = fileHashCacheActorProps,
        blacklistCache = blacklistCache)).withDispatcher(EngineDispatcher)
  }
}

/**
  * Class that orchestrates a single workflow.
  */
class WorkflowActor(workflowToStart: WorkflowToStart,
                    conf: Config,
                    callCachingEnabled: Boolean,
                    invalidateBadCacheResults:Boolean,
                    ioActor: ActorRef,
                    override val serviceRegistryActor: ActorRef,
                    workflowLogCopyRouter: ActorRef,
                    jobStoreActor: ActorRef,
                    subWorkflowStoreActor: ActorRef,
                    callCacheReadActor: ActorRef,
                    callCacheWriteActor: ActorRef,
                    dockerHashActor: ActorRef,
                    jobTokenDispenserActor: ActorRef,
                    workflowStoreActor: ActorRef,
                    backendSingletonCollection: BackendSingletonCollection,
                    serverMode: Boolean,
                    workflowHeartbeatConfig: WorkflowHeartbeatConfig,
                    totalJobsByRootWf: AtomicInteger,
                    // `Props` and not an `ActorRef` since the `RootWorkflowFileHashCacheActor` should be created as a
                    // child of this actor. The sbt subproject of `RootWorkflowFileHashCacheActor` is not visible from
                    // the subproject this class belongs to so the `Props` are passed in.
                    fileHashCacheActorProps: Option[Props],
                    blacklistCache: Option[BlacklistCache])
  extends LoggingFSM[WorkflowActorState, WorkflowActorData] with WorkflowLogging with WorkflowMetadataHelper
  with WorkflowInstrumentation with Timers {

  implicit val ec = context.dispatcher
  private val WorkflowToStart(workflowId, submissionTime, sources, initialStartableState, hogGroup) = workflowToStart
  override val workflowIdForLogging = workflowId.toPossiblyNotRoot
  override val rootWorkflowIdForLogging = workflowId.toRoot

  private val restarting = initialStartableState.restarted

  private val startTime = System.currentTimeMillis()

  private val deleteWorkflowFiles = conf.getBoolean("system.delete-workflow-files")

  private val workflowDockerLookupActor = context.actorOf(
    WorkflowDockerLookupActor.props(workflowId, dockerHashActor, initialStartableState.restarted), s"WorkflowDockerLookupActor-$workflowId")

  protected val pathBuilderFactories: List[PathBuilderFactory] = EngineFilesystems.configuredPathBuilderFactories

  startWith(WorkflowUnstartedState, WorkflowActorData(initialStartableState))

  pushCurrentStateToMetadataService(workflowId, WorkflowUnstartedState.workflowState)
  sendHeartbeat()

  // Send heartbeats possibly redundantly but unacked.
  timers.startPeriodicTimer(workflowId, SendWorkflowHeartbeatCommand, workflowHeartbeatConfig.heartbeatInterval)

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() {
    case exception if stateName == MaterializingWorkflowDescriptorState =>
      self ! MaterializeWorkflowDescriptorFailureResponse(exception)
      Stop
    case exception if stateName == InitializingWorkflowState =>
      self ! WorkflowInitializationFailedResponse(List(exception))
      Stop
    case exception if stateName == ExecutingWorkflowState =>
      self ! AbortWorkflowWithExceptionCommand(exception)
      Stop
    case exception if stateName == FinalizingWorkflowState =>
      self ! WorkflowFinalizationFailedResponse(List(exception))
      Stop
    case exception =>
      context.parent ! WorkflowFailedResponse(workflowId, stateData.lastStateReached.state, List(exception))
      context stop self
      Stop
  }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val actor = context.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryActor, workflowId, importLocalFilesystem = !serverMode, ioActorProxy = ioActor, hogGroup = hogGroup),
        "MaterializeWorkflowDescriptorActor")
      pushWorkflowStart(workflowId)
      actor ! MaterializeWorkflowDescriptorCommand(sources, conf, callCachingEnabled, invalidateBadCacheResults)
      goto(MaterializingWorkflowDescriptorState) using stateData.copy(currentLifecycleStateActor = Option(actor))
    // If the workflow is not being restarted then we can abort it immediately as nothing happened yet
    case Event(AbortWorkflowCommand, _) if !restarting => goto(WorkflowAbortedState)
  }

  /* *************************** */
  /* ****** Materializing ****** */
  /* *************************** */

  when(MaterializingWorkflowDescriptorState) {
    case Event(MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor), data) =>
      val initializerActor = context.actorOf(
        WorkflowInitializationActor.props(
          workflowIdForLogging,
          rootWorkflowIdForLogging,
          workflowDescriptor,
          ioActor,
          serviceRegistryActor,
          restarting
        ),
        name = s"WorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using data.copy(currentLifecycleStateActor = Option(initializerActor), workflowDescriptor = Option(workflowDescriptor))
    case Event(MaterializeWorkflowDescriptorFailureResponse(reason: Throwable), data) =>
      goto(WorkflowFailedState) using data.copy(lastStateReached = StateCheckpoint(MaterializingWorkflowDescriptorState, Option(List(reason))))
    // If the workflow is not being restarted then we can abort it immediately as nothing happened yet
    case Event(AbortWorkflowCommand, _) if !restarting => goto(WorkflowAbortedState)
  }

  /* ************************** */
  /* ****** Initializing ****** */
  /* ************************** */

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse(initializationData), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      val executionActor = context.actorOf(WorkflowExecutionActor.props(
        workflowDescriptor,
        ioActor = ioActor,
        serviceRegistryActor = serviceRegistryActor,
        jobStoreActor = jobStoreActor,
        subWorkflowStoreActor = subWorkflowStoreActor,
        callCacheReadActor = callCacheReadActor,
        callCacheWriteActor = callCacheWriteActor,
        workflowDockerLookupActor = workflowDockerLookupActor,
        jobTokenDispenserActor = jobTokenDispenserActor,
        backendSingletonCollection,
        initializationData,
        startState = data.effectiveStartableState,
        rootConfig = conf,
        totalJobsByRootWf = totalJobsByRootWf,
        fileHashCacheActor = fileHashCacheActorProps map context.system.actorOf,
        blacklistCache = blacklistCache), name = s"WorkflowExecutionActor-$workflowId")

      executionActor ! ExecuteWorkflowCommand
      
      val nextState = data.effectiveStartableState match {
        case RestartableAborting => WorkflowAbortingState
        case _ => ExecutingWorkflowState
      }
      goto(nextState) using data.copy(currentLifecycleStateActor = Option(executionActor), initializationData = initializationData)
    case Event(WorkflowInitializationFailedResponse(reason), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, Map.empty, CallOutputs.empty, Option(reason.toList))

    // If the workflow is not restarting, handle the Abort command normally and send an abort message to the init actor
    case Event(AbortWorkflowCommand, data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) if !restarting =>
      handleAbortCommand(data, workflowDescriptor)
  }
  
  /* ********************* */
  /* ****** Running ****** */
  /* ********************* */
  
  // Handles workflow completion events from the WEA and abort command
  val executionResponseHandler: StateFunction = {
    // Workflow responses
    case Event(WorkflowExecutionSucceededResponse(jobKeys, rootAndSubworklowIds, finalOutputs, allOutputs),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, jobKeys, finalOutputs, None, allOutputs, rootAndSubworklowIds)
    case Event(WorkflowExecutionFailedResponse(jobKeys, failures),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, jobKeys, CallOutputs.empty, Option(List(failures)))
    case Event(WorkflowExecutionAbortedResponse(jobKeys),
    data @ WorkflowActorData(_, Some(workflowDescriptor), _, StateCheckpoint(_, failures), _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, jobKeys, CallOutputs.empty, failures)

    // Whether we're running or aborting, restarting or not, pass along the abort command.
    // Note that aborting a workflow multiple times will result in as many abort commands sent to the execution actor
    case Event(AbortWorkflowWithExceptionCommand(ex), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      handleAbortCommand(data, workflowDescriptor, Option(ex))
    case Event(AbortWorkflowCommand, data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      handleAbortCommand(data, workflowDescriptor)
  }

  when(ExecutingWorkflowState)(executionResponseHandler)

  /* ********************** */
  /* ****** Aborting ****** */
  /* ********************** */

  // Handles initialization responses we can get if the abort came in when we were initializing the workflow
  val abortHandler: StateFunction = {
    // If the initialization failed, record the failure in the data and finalize the workflow
    case Event(WorkflowInitializationFailedResponse(reason), data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, Map.empty, CallOutputs.empty, Option(reason.toList))

    // Otherwise (success or abort), finalize the workflow without failures
    case Event(_: WorkflowInitializationResponse, data @ WorkflowActorData(_, Some(workflowDescriptor), _, _, _, _, _, _)) =>
      finalizeWorkflow(data, workflowDescriptor, Map.empty, CallOutputs.empty, failures = None)
  }
  
  // In aborting state, we can receive initialization responses or execution responses.
  when(WorkflowAbortingState)(abortHandler.orElse(executionResponseHandler))

  /* ************************ */
  /* ****** Finalizing ****** */
  /* ************************ */

  // When finalizing, we only expect finalization success or failure, finalization cannot be aborted.
  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, data) => finalizationSucceeded(data)
    case Event(WorkflowFinalizationFailedResponse(finalizationFailures), data) =>
      val failures = data.lastStateReached.failures.getOrElse(List.empty) ++ finalizationFailures
      goto(WorkflowFailedState) using data.copy(lastStateReached = StateCheckpoint(FinalizingWorkflowState, Option(failures)))
    case Event(AbortWorkflowCommand, _) => stay()
  }
  
  def handleAbortCommand(data: WorkflowActorData, workflowDescriptor: EngineWorkflowDescriptor, exceptionCausedAbortOpt: Option[Throwable] = None) = {
    data.currentLifecycleStateActor match {
      case Some(currentActor) =>
        currentActor ! EngineLifecycleActorAbortCommand
        goto(WorkflowAbortingState) using data.copy(lastStateReached = StateCheckpoint(stateName, exceptionCausedAbortOpt.map(List(_))))
      case None => 
        workflowLogger.warn(s"Received an abort command in state $stateName but there's no lifecycle actor associated. This is an abnormal state, finalizing the workflow anyway.")
        finalizeWorkflow(data, workflowDescriptor, Map.empty, CallOutputs.empty, None)
    }
  }

  /* **************************** */
  /* ****** Deleting files ****** */
  /* **************************** */

  // TODO: Instead of simply logging to Kibana, figure out a way to tell the user what were the errors somehow (maybe through metadata?)
  // since deletion happens only if the workflow and finalization succeeded we can directly goto Succeeded state
  when(DeletingFilesState) {
    case Event(DeleteWorkflowFilesSucceededResponse(filesNotFound, callCacheInvalidationErrors), data) =>
      workflowLogger.info(s"Successfully deleted intermediate output file(s) for root workflow $rootWorkflowIdForLogging." +
        deleteFilesAdditionalError(filesNotFound, callCacheInvalidationErrors))
      goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
    case Event(DeleteWorkflowFilesFailedResponse(errors, filesNotFound, callCacheInvalidationErrors), data) =>
      workflowLogger.info(s"Failed to delete ${errors.size} intermediate output file(s) for root workflow $rootWorkflowIdForLogging." +
        deleteFilesAdditionalError(filesNotFound, callCacheInvalidationErrors) + s" Errors: ${errors.map(ExceptionUtils.getMessage)}")
      goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
  }

  // Let these messages fall through to the whenUnhandled handler:
  when(WorkflowAbortedState) { FSM.NullFunction }
  when(WorkflowFailedState) { FSM.NullFunction }
  when(WorkflowSucceededState) { FSM.NullFunction }

  whenUnhandled {
    case Event(SendWorkflowHeartbeatCommand, _) =>
      sendHeartbeat()
      stay()
    // If the workflow is being restarted, then we have to keep going to try and reconnect to the jobs - but remember that workflow is now in abort mode
    case Event(AbortWorkflowCommand, data: WorkflowActorData) if restarting =>
      stay() using data.copy(effectiveStartableState = RestartableAborting)
    case Event(msg @ EngineStatsActor.JobCountQuery, data) =>
      data.currentLifecycleStateActor match {
        case Some(a) => a forward msg
        case None => sender ! EngineStatsActor.NoJobs // This should be impossible, but if somehow here it's technically correct
      }
      stay()
  }

  onTransition {
    case (oldState, terminalState: WorkflowActorTerminalState) =>
      // Increment counter on final transition
      setWorkflowTimePerState(terminalState.workflowState, (System.currentTimeMillis() - startTime).millis)
      workflowLogger.debug(s"transition from {} to {}. Stopping self.", arg1 = oldState, arg2 = terminalState)
      pushWorkflowEnd(workflowId)
      WorkflowProcessingEventPublishing.publish(workflowId, workflowHeartbeatConfig.cromwellId, Finished, serviceRegistryActor)
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
        def bruteForceWorkflowOptions: WorkflowOptions = sources.workflowOptions
        val system = context.system
        val ec = context.system.dispatcher
        def bruteForcePathBuilders: Future[List[PathBuilder]] = {
          // Protect against path builders that may throw an exception instead of returning a failed future
          Future(EngineFilesystems.pathBuildersForWorkflow(bruteForceWorkflowOptions, pathBuilderFactories)(system))(ec).flatten
        }

        val (workflowOptions, pathBuilders) = stateData.workflowDescriptor match {
          case Some(wd) => (wd.backendDescriptor.workflowOptions, Future.successful(wd.pathBuilders))
          case None => (bruteForceWorkflowOptions, bruteForcePathBuilders)
        }

        workflowOptions.get(FinalWorkflowLogDir).toOption match {
          case Some(destinationDir) =>
            pathBuilders
              .map(pb => workflowLogCopyRouter ! CopyWorkflowLogsActor.Copy(workflowId, PathFactory.buildPath(destinationDir, pb)))(ec)
              .recover { case e => log.error(e, "Failed to copy workflow log") }(ec)
          case None => workflowLogger.close(andDelete = WorkflowLogger.isTemporary) match {
            case Failure(f) => log.error(f, "Failed to delete workflow log")
            case _ =>
          }
        }
      }
      context stop self
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

  private def deleteFilesAdditionalError(filesNotFound: List[Path], callCacheInvalidationErrors: List[Throwable]): String = {
    val filesNotFoundMsg =
      if (filesNotFound.nonEmpty) {
        s" File(s) not found during deletion: ${filesNotFound.mkString(",")}"
      } else {
        ""
      }

    val invalidationErrorMsg =
      if (callCacheInvalidationErrors.nonEmpty) {
        s" Call cache invalidation error(s): ${callCacheInvalidationErrors.map(ExceptionUtils.getMessage)}"
      } else {
        ""
      }

    filesNotFoundMsg + invalidationErrorMsg
  }

  private def finalizationSucceeded(data: WorkflowActorData) = {
    val finalState = data.lastStateReached match {
      case StateCheckpoint(WorkflowAbortingState, None) => WorkflowAbortedState
      case StateCheckpoint(_, Some(_)) => WorkflowFailedState
      case StateCheckpoint(_, None) => WorkflowSucceededState
    }

    finalState match {
      case WorkflowSucceededState => deleteFilesOrGotoFinalState(data)
      case _ => goto(finalState) using data.copy(currentLifecycleStateActor = None)
    }
  }

  /*
    Checks if Cromwell and the user has enabled deletion of intermediate output files, and if it is enabled,
    it instantiates the DeleteWorkflowFilesActor and waits for it to respond.
    Note: We can't start deleting files before finalization succeeds as we don't want to start deleting
    them as they are being copied to another location.
  */
  private def deleteFilesOrGotoFinalState(data: WorkflowActorData) = {

    def deleteFiles() = {
      val rootWorkflowId = data.workflowDescriptor.get.rootWorkflowId
      val deleteActor = context.actorOf(DeleteWorkflowFilesActor.props(
        rootWorkflowId = rootWorkflowId,
        rootAndSubworkflowIds = data.rootAndSubworkflowIds,
        workflowFinalOutputs = data.workflowFinalOutputs,
        workflowAllOutputs = data.workflowAllOutputs,
        pathBuilders = data.workflowDescriptor.get.pathBuilders,
        serviceRegistryActor = serviceRegistryActor,
        ioActor = ioActor),
        name = s"DeleteWorkflowFilesActor-${rootWorkflowId.id}")

      deleteActor ! StartWorkflowFilesDeletion

      goto(DeletingFilesState) using data
    }

    val userDeleteFileWfOption = data.workflowDescriptor.flatMap(
      _.backendDescriptor.workflowOptions.getBoolean("delete_intermediate_output_files").toOption
    ).getOrElse(false)

    (deleteWorkflowFiles, userDeleteFileWfOption, data.workflowAllOutputs.nonEmpty) match {
      case (true, true, true) => deleteFiles()
      case (true, true, false) =>
        log.info(s"Workflow does not have any outputs. No intermediate files to delete.")
        goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
      case (true, false, _) =>
        // user has not enabled delete intermediate outputs option, so go to succeeded status
        goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
      case (false, true, _) =>
        log.info(s"User wants to delete intermediate files but it is not enabled in Cromwell config. To use it system.delete-workflow-files to true.")
        goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
      case (false, false, _) => goto(WorkflowSucceededState) using data.copy(currentLifecycleStateActor = None)
    }
  }

  private[workflow] def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs) = {
    val copyWorkflowOutputsActorProps = stateName match {
      case InitializingWorkflowState => None
      case _ => Option(CopyWorkflowOutputsActor.props(workflowIdForLogging, ioActor, workflowDescriptor, workflowOutputs, stateData.initializationData))
    }
    
    context.actorOf(WorkflowFinalizationActor.props(
      workflowDescriptor = workflowDescriptor,
      ioActor = ioActor,
      jobExecutionMap = jobExecutionMap,
      workflowOutputs = workflowOutputs,
      initializationData = stateData.initializationData,
      copyWorkflowOutputsActor = copyWorkflowOutputsActorProps
    ), name = s"WorkflowFinalizationActor")
  }

  /**
    * Run finalization actor and transition to FinalizingWorkflowState.
    */
  private def finalizeWorkflow(data: WorkflowActorData,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               jobExecutionMap: JobExecutionMap,
                               workflowFinalOutputs: CallOutputs,
                               failures: Option[List[Throwable]],
                               workflowAllOutputs: Set[WomValue] = Set.empty,
                               rootAndSubworkflowIds: Set[WorkflowId] = Set.empty) = {
    val finalizationActor = makeFinalizationActor(workflowDescriptor, jobExecutionMap, workflowFinalOutputs)
    finalizationActor ! StartFinalizationCommand
    goto(FinalizingWorkflowState) using data.copy(
      lastStateReached = StateCheckpoint (stateName, failures),
      workflowFinalOutputs = workflowFinalOutputs.outputs.values.toSet,
      workflowAllOutputs = workflowAllOutputs,
      rootAndSubworkflowIds = rootAndSubworkflowIds
    )
  }

  private def sendHeartbeat(): Unit = workflowStoreActor ! WorkflowStoreWriteHeartbeatCommand(workflowId, submissionTime)

}
