package cromwell.engine.workflow.lifecycle.execution

import java.util.concurrent.atomic.AtomicInteger

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{ActorRef, FSM, LoggingFSM, OneForOneStrategy, Props, SupervisorStrategy}
import com.typesafe.config.Config
import cromwell.backend.standard.callcaching.BlacklistCache
import cromwell.backend.{AllBackendInitializationData, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.core.logging.JobLogging
import cromwell.engine.backend.{BackendConfiguration, BackendSingletonCollection}
import cromwell.engine.workflow.WorkflowMetadataHelper
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.execution.SubWorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation.{CallPreparationFailed, Start}
import cromwell.engine.workflow.lifecycle.execution.job.preparation.SubWorkflowPreparationActor
import cromwell.engine.workflow.lifecycle.execution.job.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.engine.workflow.lifecycle.execution.keys.SubWorkflowKey
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.engine.workflow.workflowstore.StartableState
import cromwell.engine.{EngineIoFunctions, EngineWorkflowDescriptor, SubWorkflowStart}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.subworkflowstore.SubWorkflowStoreActor._
import wom.values.WomEvaluatedCallInputs

class SubWorkflowExecutionActor(key: SubWorkflowKey,
                                parentWorkflow: EngineWorkflowDescriptor,
                                expressionLanguageFunctions: EngineIoFunctions,
                                factories: Map[String, BackendLifecycleActorFactory],
                                ioActor: ActorRef,
                                override val serviceRegistryActor: ActorRef,
                                jobStoreActor: ActorRef,
                                subWorkflowStoreActor: ActorRef,
                                callCacheReadActor: ActorRef,
                                callCacheWriteActor: ActorRef,
                                workflowDockerLookupActor: ActorRef,
                                jobTokenDispenserActor: ActorRef,
                                backendSingletonCollection: BackendSingletonCollection,
                                initializationData: AllBackendInitializationData,
                                startState: StartableState,
                                rootConfig: Config,
                                totalJobsByRootWf: AtomicInteger,
                                fileHashCacheActor: Option[ActorRef],
                                blacklistCache: Option[BlacklistCache]) extends LoggingFSM[SubWorkflowExecutionActorState, SubWorkflowExecutionActorData] with JobLogging with WorkflowMetadataHelper with CallMetadataHelper {

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  override val workflowIdForLogging = parentWorkflow.possiblyNotRootWorkflowId
  override val rootWorkflowIdForLogging = parentWorkflow.rootWorkflowId
  override val workflowIdForCallMetadata = parentWorkflow.id
  override def jobTag: String = key.tag

  startWith(SubWorkflowPendingState, SubWorkflowExecutionActorData.empty)

  override def preStart(): Unit = {
    context.system.eventStream.publish(SubWorkflowStart(self))
    super.preStart()
  }

  private var eventList: Seq[ExecutionEvent] = Seq(ExecutionEvent(stateName.toString))

  when(SubWorkflowPendingState) {
    case Event(Execute, _) =>
      if (startState.restarted) {
        subWorkflowStoreActor ! QuerySubWorkflow(parentWorkflow.id, key)
        goto(SubWorkflowCheckingStoreState)
      } else {
        requestValueStore(createSubWorkflowId())
      }
  }

  when(SubWorkflowCheckingStoreState) {
    case Event(SubWorkflowFound(entry), _) =>
      requestValueStore(WorkflowId.fromString(entry.subWorkflowExecutionUuid))
    case Event(_: SubWorkflowNotFound, _) =>
      requestValueStore(createSubWorkflowId())
    case Event(SubWorkflowStoreFailure(command, reason), _) =>
      jobLogger.error(reason, s"SubWorkflowStore failure for command $command, starting sub workflow with fresh ID.")
      requestValueStore(createSubWorkflowId())
  }

  /*
    * ! Hot Potato Warning !
    * We ask explicitly for the output store so we can use it on the fly and more importantly not store it as a
    * variable in this actor, which would prevent it from being garbage collected for the duration of the
    * subworkflow and would lead to memory leaks.
    */
  when(WaitingForValueStore) {
    case Event(valueStore: ValueStore, SubWorkflowExecutionActorData(Some(subWorkflowId), _)) =>
      prepareSubWorkflow(subWorkflowId, valueStore)
    case Event(_: ValueStore, _) =>
      context.parent ! SubWorkflowFailedResponse(key, Map.empty, new IllegalStateException(
        "This is a programmer error, we're ready to prepare the job and should have" +
          " a SubWorkflowId to use by now but somehow haven't. Failing the workflow."))
      context stop self
      stay()
  }

  when(SubWorkflowPreparingState) {
    case Event(SubWorkflowPreparationSucceeded(subWorkflowEngineDescriptor, inputs), data) =>
      startSubWorkflow(subWorkflowEngineDescriptor, inputs, data)
    case Event(failure: CallPreparationFailed, _) =>
      context.parent ! SubWorkflowFailedResponse(key, Map.empty, failure.throwable)
      context stop self
      stay()
  }

  when(SubWorkflowRunningState) {
    case Event(WorkflowExecutionSucceededResponse(executedJobKeys, outputs), _) =>
      context.parent ! SubWorkflowSucceededResponse(key, executedJobKeys, outputs)
      goto(SubWorkflowSucceededState)
    case Event(WorkflowExecutionFailedResponse(executedJobKeys, reason), _) =>
      context.parent ! SubWorkflowFailedResponse(key, executedJobKeys, reason)
      goto(SubWorkflowFailedState)
    case Event(WorkflowExecutionAbortedResponse(executedJobKeys), _) =>
      context.parent ! SubWorkflowAbortedResponse(key, executedJobKeys)
      goto(SubWorkflowAbortedState)
    case Event(EngineLifecycleActorAbortCommand, SubWorkflowExecutionActorData(_, Some(actorRef))) =>
      actorRef ! EngineLifecycleActorAbortCommand
      stay()
  }

  when(SubWorkflowSucceededState) { FSM.NullFunction }
  when(SubWorkflowFailedState) { FSM.NullFunction }
  when(SubWorkflowAbortedState) { FSM.NullFunction }

  whenUnhandled {
    case Event(SubWorkflowStoreRegisterSuccess(_), _) =>
      // Nothing to do here
      stay()
    case Event(SubWorkflowStoreFailure(command, reason), _) =>
      jobLogger.error(reason, s"SubWorkflowStore failure for command $command")
      stay()
    case Event(EngineLifecycleActorAbortCommand, _) =>
      context.parent ! SubWorkflowAbortedResponse(key, Map.empty)
      goto(SubWorkflowAbortedState)
  }

  onTransition {
    case (_, toState) =>
      stateData.subWorkflowId foreach { id => pushCurrentStateToMetadataService(id, toState.workflowState) }
  }

  onTransition {
    case (_, _: SubWorkflowTerminalState) =>
      stateData.subWorkflowId match {
        case Some(id) =>
          pushWorkflowEnd(id)
          pushExecutionEventsToMetadataService(key, eventList)
        case None => jobLogger.error("Sub workflow completed without a Sub Workflow UUID.")
      }
      context stop self
  }

  onTransition {
    case _ -> toState => eventList :+= ExecutionEvent(toState.toString)
  }

  private def startSubWorkflow(subWorkflowEngineDescriptor: EngineWorkflowDescriptor, inputs: WomEvaluatedCallInputs, data: SubWorkflowExecutionActorData) = {
    val subWorkflowActor = createSubWorkflowActor(subWorkflowEngineDescriptor)

    subWorkflowActor ! WorkflowExecutionActor.ExecuteWorkflowCommand
    context.parent ! JobRunning(key, inputs)
    pushWorkflowRunningMetadata(subWorkflowEngineDescriptor.backendDescriptor, inputs)

    goto(SubWorkflowRunningState) using data.copy(subWorkflowActor = Option(subWorkflowActor))
  }

  private def prepareSubWorkflow(subWorkflowId: WorkflowId, valueStore: ValueStore) = {
    createSubWorkflowPreparationActor(subWorkflowId) ! Start(valueStore)
    context.parent ! JobStarting(key)
    pushCurrentStateToMetadataService(subWorkflowId, WorkflowRunning)
    pushWorkflowStart(subWorkflowId)
    goto(SubWorkflowPreparingState) using SubWorkflowExecutionActorData(Option(subWorkflowId), None)
  }

  private def requestValueStore(workflowId: WorkflowId) = {
    context.parent ! RequestValueStore
    goto(WaitingForValueStore) using SubWorkflowExecutionActorData(Option(workflowId), None)
  }

  def createSubWorkflowPreparationActor(subWorkflowId: WorkflowId) = {
    context.actorOf(
      SubWorkflowPreparationActor.props(parentWorkflow, expressionLanguageFunctions, key, subWorkflowId),
      s"$subWorkflowId-SubWorkflowPreparationActor-${key.tag}"
    )
  }

  def createSubWorkflowActor(subWorkflowEngineDescriptor: EngineWorkflowDescriptor) = {
    context.actorOf(
      WorkflowExecutionActor.props(
        subWorkflowEngineDescriptor,
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
        startState,
        rootConfig,
        totalJobsByRootWf,
        fileHashCacheActor = fileHashCacheActor,
        blacklistCache = blacklistCache
      ),
      s"${subWorkflowEngineDescriptor.id}-SubWorkflowActor-${key.tag}"
    )
  }

  private def pushWorkflowRunningMetadata(subWorkflowDescriptor: BackendWorkflowDescriptor, workflowInputs: WomEvaluatedCallInputs) = {
    val subWorkflowId = subWorkflowDescriptor.id
    val parentWorkflowMetadataKey = MetadataKey(parentWorkflow.id, Option(MetadataJobKey(key.node.fullyQualifiedName, key.index, key.attempt)), CallMetadataKeys.SubWorkflowId)

    val events = List(
      MetadataEvent(parentWorkflowMetadataKey, MetadataValue(subWorkflowId)),
      MetadataEvent(MetadataKey(subWorkflowId, None, WorkflowMetadataKeys.Name), MetadataValue(key.node.localName)),
      MetadataEvent(MetadataKey(subWorkflowId, None, WorkflowMetadataKeys.ParentWorkflowId), MetadataValue(parentWorkflow.id)),
      MetadataEvent(MetadataKey(subWorkflowId, None, WorkflowMetadataKeys.RootWorkflowId), MetadataValue(parentWorkflow.rootWorkflow.id))
    )

    val inputEvents = workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(subWorkflowId, None,WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (inputName, womValue) =>
          womValueToMetadataEvents(MetadataKey(subWorkflowId, None, s"${WorkflowMetadataKeys.Inputs}:${inputName.name}"), womValue)
        }
    }

    val workflowRootEvents = buildWorkflowRootMetadataEvents(subWorkflowDescriptor)

    serviceRegistryActor ! PutMetadataAction(events ++ inputEvents ++ workflowRootEvents)
  }

  private def buildWorkflowRootMetadataEvents(subWorkflowDescriptor: BackendWorkflowDescriptor) = {
    val subWorkflowId = subWorkflowDescriptor.id

    factories flatMap {
      case (backendName, factory) =>
        BackendConfiguration.backendConfigurationDescriptor(backendName).toOption map { config =>
          backendName -> factory.getWorkflowExecutionRootPath(subWorkflowDescriptor, config.backendConfig, initializationData.get(backendName))
        }
    } map {
      case (backend, wfRoot) =>
        MetadataEvent(MetadataKey(subWorkflowId, None, s"${WorkflowMetadataKeys.WorkflowRoot}[$backend]"), MetadataValue(wfRoot.toAbsolutePath))
    }
  }

  private def createSubWorkflowId() = {
    val subWorkflowId = WorkflowId.randomId()
    // Register ID to the sub workflow store
    subWorkflowStoreActor ! RegisterSubWorkflow(parentWorkflow.rootWorkflow.id, parentWorkflow.id, key, subWorkflowId)
    subWorkflowId
  }
}

object SubWorkflowExecutionActor {
  sealed trait SubWorkflowExecutionActorState {
    def workflowState: WorkflowState
  }
  sealed trait SubWorkflowTerminalState extends SubWorkflowExecutionActorState

  case object SubWorkflowPendingState extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowRunning
  }
  case object SubWorkflowCheckingStoreState extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowRunning
  }
  case object SubWorkflowPreparingState extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowRunning
  }
  case object SubWorkflowRunningState extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowRunning
  }
  case object WaitingForValueStore extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowRunning
  }
  case object SubWorkflowAbortingState extends SubWorkflowExecutionActorState {
    override val workflowState = WorkflowAborting
  }

  case object SubWorkflowSucceededState extends SubWorkflowTerminalState {
    override val workflowState = WorkflowSucceeded
  }
  case object SubWorkflowAbortedState extends SubWorkflowTerminalState {
    override val workflowState = WorkflowAborted
  }
  case object SubWorkflowFailedState extends SubWorkflowTerminalState {
    override val workflowState = WorkflowFailed
  }

  object SubWorkflowExecutionActorData {
    def empty = SubWorkflowExecutionActorData(None, None)
  }
  case class SubWorkflowExecutionActorData(subWorkflowId: Option[WorkflowId], subWorkflowActor: Option[ActorRef])

  sealed trait EngineWorkflowExecutionActorCommand
  case object Execute

  def props(key: SubWorkflowKey,
            parentWorkflow: EngineWorkflowDescriptor,
            expressionLanguageFunctions: EngineIoFunctions,
            factories: Map[String, BackendLifecycleActorFactory],
            ioActor: ActorRef,
            serviceRegistryActor: ActorRef,
            jobStoreActor: ActorRef,
            subWorkflowStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            callCacheWriteActor: ActorRef,
            workflowDockerLookupActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection,
            initializationData: AllBackendInitializationData,
            startState: StartableState,
            rootConfig: Config,
            totalJobsByRootWf: AtomicInteger,
            fileHashCacheActor: Option[ActorRef],
            blacklistCache: Option[BlacklistCache]) = {
    Props(new SubWorkflowExecutionActor(
      key,
      parentWorkflow,
      expressionLanguageFunctions,
      factories,
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
      startState,
      rootConfig,
      totalJobsByRootWf,
      fileHashCacheActor = fileHashCacheActor,
      blacklistCache = blacklistCache)
    ).withDispatcher(EngineDispatcher)
  }
}
