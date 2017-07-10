package cromwell.engine.workflow.lifecycle.execution

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{ActorRef, FSM, LoggingFSM, OneForOneStrategy, Props, SupervisorStrategy}
import cromwell.backend.{AllBackendInitializationData, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import cromwell.core._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.logging.JobLogging
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, BackendSingletonCollection}
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{CallPreparationFailed, Start}
import cromwell.engine.workflow.lifecycle.execution.SubWorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.preparation.SubWorkflowPreparationActor
import cromwell.engine.workflow.lifecycle.execution.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.subworkflowstore.SubWorkflowStoreActor._
import wdl4s.EvaluatedTaskInputs

class SubWorkflowExecutionActor(key: SubWorkflowKey,
                                data: WorkflowExecutionActorData,
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
                                restarting: Boolean) extends LoggingFSM[SubWorkflowExecutionActorState, SubWorkflowExecutionActorData] with JobLogging with WorkflowMetadataHelper with CallMetadataHelper {

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  private val parentWorkflow = data.workflowDescriptor
  override val workflowId = parentWorkflow.id
  override val workflowIdForCallMetadata = parentWorkflow.id
  override def jobTag: String = key.tag

  startWith(SubWorkflowPendingState, SubWorkflowExecutionActorData.empty)

  private var eventList: Seq[ExecutionEvent] = Seq(ExecutionEvent(stateName.toString))

  when(SubWorkflowPendingState) {
    case Event(Execute, _) =>
      if (restarting) {
        subWorkflowStoreActor ! QuerySubWorkflow(parentWorkflow.id, key)
        goto(SubWorkflowCheckingStoreState)
      } else {
        prepareSubWorkflow(createSubWorkflowId())
      }
  }

  when(SubWorkflowCheckingStoreState) {
    case Event(SubWorkflowFound(entry), _) =>
      prepareSubWorkflow(WorkflowId.fromString(entry.subWorkflowExecutionUuid))
    case Event(_: SubWorkflowNotFound, _) =>
      prepareSubWorkflow(createSubWorkflowId())
    case Event(SubWorkflowStoreFailure(command, reason), _) =>
      jobLogger.error(reason, s"SubWorkflowStore failure for command $command, starting sub workflow with fresh ID.")
      prepareSubWorkflow(createSubWorkflowId())
  }

  when(SubWorkflowPreparingState) {
    case Event(SubWorkflowPreparationSucceeded(subWorkflowEngineDescriptor, inputs), _) =>
      startSubWorkflow(subWorkflowEngineDescriptor, inputs)
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

  private def startSubWorkflow(subWorkflowEngineDescriptor: EngineWorkflowDescriptor, inputs: EvaluatedTaskInputs) = {
    val subWorkflowActor = createSubWorkflowActor(subWorkflowEngineDescriptor)

    subWorkflowActor ! WorkflowExecutionActor.ExecuteWorkflowCommand
    context.parent ! JobRunning(key, inputs, Option(subWorkflowActor))
    pushWorkflowRunningMetadata(subWorkflowEngineDescriptor.backendDescriptor, inputs)

    goto(SubWorkflowRunningState)
  }

  private def prepareSubWorkflow(subWorkflowId: WorkflowId) = {
    createSubWorkflowPreparationActor(subWorkflowId) ! Start
    context.parent ! JobStarting(key)
    pushCurrentStateToMetadataService(subWorkflowId, WorkflowRunning)
    pushWorkflowStart(subWorkflowId)
    goto(SubWorkflowPreparingState) using SubWorkflowExecutionActorData(Option(subWorkflowId))
  }

  def createSubWorkflowPreparationActor(subWorkflowId: WorkflowId) = {
    context.actorOf(
      SubWorkflowPreparationActor.props(data, key, subWorkflowId),
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
      restarting
    ),
      s"${subWorkflowEngineDescriptor.id}-SubWorkflowActor-${key.tag}"
    )
  }

  private def pushWorkflowRunningMetadata(subWorkflowDescriptor: BackendWorkflowDescriptor, workflowInputs: EvaluatedTaskInputs) = {
    val subWorkflowId = subWorkflowDescriptor.id
    val parentWorkflowMetadataKey = MetadataKey(parentWorkflow.id, Option(MetadataJobKey(key.scope.fullyQualifiedName, key.index, key.attempt)), CallMetadataKeys.SubWorkflowId)
    
    val events = List(
      MetadataEvent(parentWorkflowMetadataKey, MetadataValue(subWorkflowId)),
      MetadataEvent(MetadataKey(subWorkflowId, None, WorkflowMetadataKeys.Name), MetadataValue(key.scope.callable.unqualifiedName)),
      MetadataEvent(MetadataKey(subWorkflowId, None, WorkflowMetadataKeys.ParentWorkflowId), MetadataValue(parentWorkflow.id))
    )

    val inputEvents = workflowInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(MetadataKey(subWorkflowId, None,WorkflowMetadataKeys.Inputs)))
      case inputs =>
        inputs flatMap { case (inputName, wdlValue) =>
          wdlValueToMetadataEvents(MetadataKey(subWorkflowId, None, s"${WorkflowMetadataKeys.Inputs}:${inputName.unqualifiedName}"), wdlValue)
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
    def empty = SubWorkflowExecutionActorData(None)
  }
  case class SubWorkflowExecutionActorData(subWorkflowId: Option[WorkflowId])

  sealed trait EngineWorkflowExecutionActorCommand
  case object Execute

  def props(key: SubWorkflowKey,
            data: WorkflowExecutionActorData,
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
            restarting: Boolean) = {
    Props(new SubWorkflowExecutionActor(
      key,
      data,
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
      restarting)
    ).withDispatcher(EngineDispatcher)
  }
}