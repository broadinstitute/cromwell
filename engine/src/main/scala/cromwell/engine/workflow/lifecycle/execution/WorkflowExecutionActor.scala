package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.ExecutionIndex._
import cromwell.core.ExecutionStatus._
import cromwell.core.ExecutionStore.ExecutionStoreEntry
import cromwell.core.OutputStore.OutputEntry
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{WorkflowId, _}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.CromwellBackends
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.WorkflowExecutionActorState
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}
import cromwell.services.MetadataServiceActor._
import cromwell.services._
import lenthall.exception.ThrowableAggregation
import wdl4s._
import wdl4s.types.WdlArrayType
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlValue}

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.util.{Failure, Random, Success, Try}

object WorkflowExecutionActor {

  /**
    * States
    */
  sealed trait WorkflowExecutionActorState { def terminal = false }
  sealed trait WorkflowExecutionActorTerminalState extends WorkflowExecutionActorState { override val terminal = true }

  case object WorkflowExecutionPendingState extends WorkflowExecutionActorState
  case object WorkflowExecutionInProgressState extends WorkflowExecutionActorState
  case object WorkflowExecutionAbortingState extends WorkflowExecutionActorState
  case object WorkflowExecutionSuccessfulState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionFailedState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionAbortedState extends WorkflowExecutionActorTerminalState

  /**
    * Commands
    */
  sealed trait WorkflowExecutionActorCommand
  case object StartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object RestartExecutingWorkflowCommand extends WorkflowExecutionActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowExecutionActorResponse {
    def executionStore: ExecutionStore

    def outputStore: OutputStore
  }

  case class WorkflowExecutionSucceededResponse(executionStore: ExecutionStore, outputStore: OutputStore)
    extends WorkflowExecutionActorResponse {
    override def toString = "WorkflowExecutionSucceededResponse"
  }

  case class WorkflowExecutionAbortedResponse(executionStore: ExecutionStore, outputStore: OutputStore)
    extends WorkflowExecutionActorResponse with EngineLifecycleActorAbortedResponse {
    override def toString = "WorkflowExecutionAbortedResponse"
  }

  final case class WorkflowExecutionFailedResponse(executionStore: ExecutionStore, outputStore: OutputStore,
                                                   reasons: Seq[Throwable]) extends WorkflowExecutionActorResponse {
    override def toString = "WorkflowExecutionFailedResponse"
  }

  /**
    * Internal control flow messages
    */
  private case class JobInitializationFailed(jobKey: JobKey, throwable: Throwable)
  private case class ScatterCollectionFailedResponse(collectorKey: CollectorKey, throwable: Throwable)
  private case class ScatterCollectionSucceededResponse(collectorKey: CollectorKey, outputs: JobOutputs)

  /**
    * Internal ADTs
    */
  case class ScatterKey(scope: Scatter) extends JobKey {
    override val index = None // When scatters are nested, this might become Some(_)
    override val attempt = 1
    override val tag = scope.unqualifiedName

    /**
      * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
      *
      * @param count Number of ways to scatter the children.
      * @return ExecutionStore of scattered children.
      */
    def populate(count: Int): Map[JobKey, ExecutionStatus.Value] = {
      val keys = this.scope.children flatMap { explode(_, count) }
      keys map { _ -> ExecutionStatus.NotStarted } toMap
    }

    private def explode(scope: Scope, count: Int): Seq[JobKey] = {
      scope match {
        case call: Call =>
          val shards = (0 until count) map { i => BackendJobDescriptorKey(call, Option(i), 1) }
          shards :+ CollectorKey(call)
        case scatter: Scatter =>
          throw new UnsupportedOperationException("Nested Scatters are not supported (yet).")
        case e =>
          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
      }
    }
  }

  // Represents a scatter collection for a call in the execution store
  case class CollectorKey(scope: Call) extends JobKey {
    override val index = None
    override val attempt = 1
    override val tag = s"Collector-${scope.unqualifiedName}"
  }

  case class WorkflowExecutionException(override val throwables: List[Throwable]) extends ThrowableAggregation {
    override val exceptionContext = s"WorkflowExecutionActor"
  }

  def props(workflowId: WorkflowId,
            workflowDescriptor: EngineWorkflowDescriptor,
            serviceRegistryActor: ActorRef): Props = Props(WorkflowExecutionActor(workflowId, workflowDescriptor, serviceRegistryActor))

  private implicit class EnhancedExecutionStore(val executionStore: ExecutionStore) extends AnyVal {
    /** Find currently runnable scopes */
    def runnableScopes: Iterable[JobKey] = executionStore.store.filter(isRunnable).keys

    private def isRunnable(entry: ExecutionStoreEntry) = {
      entry match {
        case (key, ExecutionStatus.NotStarted) => arePrerequisitesDone(key)
        case _ => false
      }
    }

    def findShardEntries(key: CollectorKey): Iterable[ExecutionStoreEntry] = executionStore.store collect {
      case (k: BackendJobDescriptorKey, v) if k.scope == key.scope && k.isShard => (k, v)
    }

    private def arePrerequisitesDone(key: JobKey): Boolean = {
      def isDone(e: JobKey): Boolean = executionStore.store exists {
        case (k, s) => k.scope == e.scope && k.index == e.index && s == ExecutionStatus.Done
      }

      val upstream = key.scope.prerequisiteScopes.map(s => upstreamEntries(key, s))
      val downstream = key match {
        case collector: CollectorKey => findShardEntries(collector)
        case _ => Nil
      }
      val dependencies = upstream.flatten ++ downstream
      val dependenciesResolved = executionStore.store.filter(dependencies).keys forall isDone

      /**
        * We need to make sure that all prerequisiteScopes have been resolved to some entry before going forward.
        * If a scope cannot be resolved it may be because it is in a scatter that has not been populated yet,
        * therefore there is no entry in the executionStore for this scope.
        * If that's the case this prerequisiteScope has not been run yet, hence the (upstream forall {_.nonEmpty})
        */
      (upstream forall { _.nonEmpty }) && dependenciesResolved
    }

    private def upstreamEntries(entry: JobKey, prerequisiteScope: Scope): Seq[ExecutionStoreEntry] = {
      prerequisiteScope.closestCommonAncestor(entry.scope) match {
        /**
          * If this entry refers to a Scope which has a common ancestor with prerequisiteScope
          * and that common ancestor is a Scatter block, then find the shard with the same index
          * as 'entry'.  In other words, if you're in the same scatter block as your pre-requisite
          * scope, then depend on the shard (with same index).
          *
          * NOTE: this algorithm was designed for ONE-LEVEL of scattering and probably does not
          * work as-is for nested scatter blocks
          */
        case Some(ancestor: Scatter) =>
          executionStore.store filter {
            case (k, _) => k.scope == prerequisiteScope && k.index == entry.index
          } toSeq

        /**
          * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
          * on every shard of the pre-requisite scope to finish.
          */
        case _ =>
          executionStore.store filter {
            case (k, _) => k.scope == prerequisiteScope && k.index.isEmpty
          } toSeq
      }
    }
  }

  private implicit class EnhancedOutputStore(val outputStore: OutputStore) extends AnyVal {
    /**
      * Try to generate output for a collector call, by collecting outputs for all of its shards.
      * It's fail-fast on shard output retrieval
      */
    def generateCollectorOutput(collector: CollectorKey,
                                shards: Iterable[BackendJobDescriptorKey]): Try[JobOutputs] = Try {
      val shardsOutputs = shards.toSeq sortBy { _.index.fromIndex } map { e =>
        outputStore.fetchCallOutputEntries(e.scope, e.index) map {
          _.outputs
        } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
      }
      collector.scope.task.outputs map { taskOutput =>
        val wdlValues = shardsOutputs.map(
          _.getOrElse(taskOutput.name, throw new RuntimeException(s"Could not retrieve output ${taskOutput.name}")))
        val arrayOfValues = new WdlArray(WdlArrayType(taskOutput.wdlType), wdlValues)
        taskOutput.name -> JobOutput(arrayOfValues, None)
      } toMap
    }
  }

}

final case class WorkflowExecutionActor(workflowId: WorkflowId,
                                        workflowDescriptor: EngineWorkflowDescriptor,
                                        serviceRegistryActor: ActorRef)
  extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] with WorkflowLogging {

  import WorkflowExecutionActor._
  import lenthall.config.ScalaConfig._

  val tag = s"WorkflowExecutionActor [UUID(${workflowId.shortString})]"
  private lazy val DefaultMaxRetriesFallbackValue = 10

  implicit val actorSystem = context.system
  implicit val ec = context.dispatcher

  val MaxRetries = ConfigFactory.load().getIntOption("system.max-retries") match {
    case Some(value) => value
    case None =>
      workflowLogger.warn(s"Failed to load the max-retries value from the configuration. Defaulting back to a value of '$DefaultMaxRetriesFallbackValue'.")
      DefaultMaxRetriesFallbackValue
  }

  private val factories = TryUtil.sequenceMap(workflowDescriptor.backendAssignments.values.toSet[String] map { backendName =>
    backendName -> CromwellBackends.shadowBackendLifecycleFactory(backendName)
  } toMap) recover {
    case e => throw new RuntimeException("Could not instantiate backend factories", e)
  } get

  // Initialize the StateData with ExecutionStore (all calls as NotStarted) and SymbolStore
  startWith(
    WorkflowExecutionPendingState,
    WorkflowExecutionActorData(
      workflowDescriptor,
      executionStore = buildInitialExecutionStore(),
      backendJobExecutionActors = Map.empty,
      outputStore = OutputStore.empty
    )
  )

  private def buildInitialExecutionStore(): ExecutionStore = {
    val workflow = workflowDescriptor.backendDescriptor.workflowNamespace.workflow
    // Only add direct children to the store, the rest is dynamically created when necessary
    val keys = workflow.children map {
      case call: Call => BackendJobDescriptorKey(call, None, 1)
      case scatter: Scatter => ScatterKey(scatter)
    }

    ExecutionStore(keys.map(_ -> NotStarted).toMap)
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"${jobDescriptor.descriptor.id}-BackendExecutionActor-${jobDescriptor.key.tag}"
  }

  when(WorkflowExecutionPendingState) {
    case Event(StartExecutingWorkflowCommand, stateData) =>
      val data = startRunnableScopes(stateData)
      goto(WorkflowExecutionInProgressState) using data
    case Event(RestartExecutingWorkflowCommand, _) =>
      // TODO: Restart executing
      goto(WorkflowExecutionInProgressState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, actorProps), stateData) =>
      pushPreparedJobMetadata(jobDescriptor.key, jobDescriptor.inputs)
      val backendExecutionActor = context.actorOf(actorProps, buildJobExecutionActorName(jobDescriptor))
      backendExecutionActor ! ExecuteJobCommand
      pushRunningJobMetadata(jobDescriptor.key)
      stay() using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobDescriptor.key -> ExecutionStatus.Running)))
        .addBackendJobExecutionActor(jobDescriptor.key, backendExecutionActor)
    case Event(BackendJobPreparationFailed(jobKey, throwable), stateData) =>
      workflowLogger.error(throwable, "Failed to start job {}", jobKey) // TODO: This log is a candidate for removal. It's now recorded in metadata
      pushFailedJobMetadata(jobKey, None, throwable, retryableFailure = false)
      context.parent ! WorkflowExecutionFailedResponse(stateData.executionStore, stateData.outputStore, List(throwable))
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
    case Event(SucceededResponse(jobKey, returnCode, callOutputs), stateData) =>
      pushSuccessfulJobMetadata(jobKey, returnCode, callOutputs)
      handleJobSuccessful(jobKey, callOutputs, stateData)
    case Event(FailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      workflowLogger.warn(s"Job ${jobKey.tag} failed! Reason: ${reason.getMessage}") // TODO: PBE This log is a candidate for removal. It's now recorded in metadata
      pushFailedJobMetadata(jobKey, returnCode, reason, retryableFailure = false)
      context.parent ! WorkflowExecutionFailedResponse(stateData.executionStore, stateData.outputStore, List(reason))
      val mergedStateData = stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
      goto(WorkflowExecutionFailedState) using mergedStateData.removeBackendJobExecutionActor(jobKey)
    case Event(FailedRetryableResponse(jobKey, reason, returnCode), stateData) =>
      workflowLogger.warn(s"Job ${jobKey.tag} failed with a retryable failure: ${reason.getMessage}")
      pushFailedJobMetadata(jobKey, None, reason, retryableFailure = true)
      handleRetryableFailure(jobKey, reason, returnCode)
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      workflowLogger.warn(s"Jobs failed to initialize: $reason") // TODO: PBE This log is a candidate for removal. It's now recorded in metadata
      pushFailedJobMetadata(jobKey, None, reason, retryableFailure = false)
      goto(WorkflowExecutionFailedState)
    case Event(ScatterCollectionSucceededResponse(jobKey, callOutputs), stateData) =>
      handleJobSuccessful(jobKey, callOutputs, stateData)
  }

  when(WorkflowExecutionSuccessfulState) {
    FSM.NullFunction
  }
  when(WorkflowExecutionFailedState) {
    alreadyFailedMopUp
  }
  when(WorkflowExecutionAbortedState) {
    alreadyFailedMopUp
  }

  /**
    * Mop up function to handle a set of incoming results if this workflow has already failed:
    */
  private def alreadyFailedMopUp: StateFunction = {
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      pushFailedJobMetadata(jobKey, None, reason, retryableFailure = false)
      stay
    case Event(FailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedJobMetadata(jobKey, returnCode, reason, retryableFailure = false)
      stay
    case Event(FailedRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedJobMetadata(jobKey, returnCode, reason, retryableFailure = true)
      stay
    case Event(SucceededResponse(jobKey, returnCode, callOutputs), stateData) =>
      pushSuccessfulJobMetadata(jobKey, returnCode, callOutputs)
      stay
  }
  when(WorkflowExecutionAbortingState) {
    case Event(AbortedResponse(jobKey), stateData) =>
      workflowLogger.info(s"$tag job aborted: ${jobKey.tag}")
      val newStateData = stateData.removeBackendJobExecutionActor(jobKey)
      if (newStateData.backendJobExecutionActors.isEmpty) {
        workflowLogger.info(s"$tag all jobs aborted")
        goto(WorkflowExecutionAbortedState)
      } else {
        stay() using newStateData
      }
  }

  whenUnhandled {
    case Event(MetadataPutFailed(action, error), _) =>
      // Do something useful here??
      workflowLogger.warn(s"$tag Put failed for Metadata action $action : ${error.getMessage}")
      stay
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      if (stateData.backendJobExecutionActors.nonEmpty) {
        stateData.backendJobExecutionActors.values foreach {_ ! AbortJobCommand}
        goto(WorkflowExecutionAbortingState)
      } else {
        goto(WorkflowExecutionAbortedState)
      }
    case unhandledMessage =>
      workflowLogger.warn(s"$tag received an unhandled message: ${unhandledMessage.event} in state: $stateName")
      stay
  }

  onTransition {
    case fromState -> toState if toState.terminal =>
      workflowLogger.info(s"$tag transitioning from $fromState to $toState. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      workflowLogger.info(s"$tag transitioning from $fromState to $toState.")
  }

  onTransition {
    case _ -> WorkflowExecutionSuccessfulState =>
      pushWorkflowOutputMetadata(nextStateData)
      context.parent ! WorkflowExecutionSucceededResponse(nextStateData.executionStore, nextStateData.outputStore)
    case _ -> WorkflowExecutionAbortedState =>
      context.parent ! WorkflowExecutionAbortedResponse(nextStateData.executionStore, nextStateData.outputStore)
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      workflowLogger.info(s"$tag done. Shutting down.")
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey, reason: Throwable, returnCode: Option[Int]) = {
    // We start with index 1 for #attempts, hence invariant breaks only if jobKey.attempt > MaxRetries
    if (jobKey.attempt <= MaxRetries) {
      val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
      workflowLogger.info(s"Retrying job execution for ${newJobKey.tag}")
      /** Currently, we update the status of the old key to Preempted, and add a new entry (with the #attempts incremented by 1)
        * to the execution store with status as NotStarted. This allows startRunnableCalls to re-execute this job */
      val executionDiff = WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Preempted, newJobKey -> ExecutionStatus.NotStarted))
      val newData = stateData.mergeExecutionDiff(executionDiff)
      stay() using startRunnableScopes(newData)
    } else {
      workflowLogger.warn(s"Exhausted maximum number of retries for job ${jobKey.tag}. Failing.")
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
    }
  }

  private def handleJobSuccessful(jobKey: JobKey, outputs: JobOutputs, data: WorkflowExecutionActorData) = {
    workflowLogger.info(s"Job ${jobKey.tag} succeeded! Outputs: ${outputs.mkString("\n")}")
    val newData = data.jobExecutionSuccess(jobKey, outputs)

    if (newData.isWorkflowComplete) {
      workflowLogger.info(newData.outputsJson())
      goto(WorkflowExecutionSuccessfulState) using newData
    }
    else
      stay() using startRunnableScopes(newData)
  }

  private def pushWorkflowOutputMetadata(data: WorkflowExecutionActorData) = {
    val reportableOutputs = workflowDescriptor.backendDescriptor.workflowNamespace.workflow.outputs
    val keyValues = data.outputStore.store filterKeys {
      _.index.isEmpty
    } flatMap {
      case (key, value) =>
        value collect {
          case entry if isReportableOutput(key.call, entry, reportableOutputs) =>
            s"${key.call.fullyQualifiedName}.${entry.name}" -> entry.wdlValue
        }
    } collect {
      case (key, Some(wdlValue)) => (key, wdlValue)
    }

    val events = keyValues match {
      case empty if empty.isEmpty => List(MetadataEvent.empty(MetadataKey(workflowId, None, WorkflowMetadataKeys.Outputs)))
      case _ => keyValues flatMap {
        case (outputName, outputValue) =>
          wdlValueToMetadataEvents(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Outputs}:$outputName"), outputValue)
      }
    }

    serviceRegistryActor ! PutMetadataAction(events)
  }

  private def isReportableOutput(scope: Scope, entry: OutputEntry,
                                 reportableOutputs: Seq[ReportableSymbol]): Boolean = {
    reportableOutputs exists { reportableOutput =>
      reportableOutput.fullyQualifiedName == s"${scope.fullyQualifiedName}.${entry.name}"
    }
  }

  private def pushRunningJobMetadata(jobKey: JobKey) = {
    serviceRegistryActor ! PutMetadataAction(MetadataEvent(metadataKey(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Running)))
  }

  private def pushSuccessfulJobMetadata(jobKey: JobKey, returnCode: Option[Int], outputs: JobOutputs) = {
    val completionEvents = completedJobMetadataEvents(jobKey, ExecutionStatus.Done, returnCode)

    val outputEvents = outputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKey(jobKey, s"${CallMetadataKeys.Outputs}")))
      case _ =>
        outputs flatMap { case (lqn, value) => wdlValueToMetadataEvents(metadataKey(jobKey, s"${CallMetadataKeys.Outputs}:$lqn"), value.wdlValue) }
    }

    serviceRegistryActor ! PutMetadataAction(completionEvents ++ outputEvents)
  }

  private def pushFailedJobMetadata(jobKey: JobKey, returnCode: Option[Int], failure: Throwable, retryableFailure: Boolean) = {
    val completionEvents = completedJobMetadataEvents(jobKey, ExecutionStatus.Failed, returnCode)
    val retryableFailureEvent = MetadataEvent(metadataKey(jobKey, CallMetadataKeys.RetryableFailure), MetadataValue(retryableFailure))
    val failureEvents = throwableToMetadataEvents(metadataKey(jobKey, s"${CallMetadataKeys.Failures}[$randomNumberString]"), failure).+:(retryableFailureEvent)

    serviceRegistryActor ! PutMetadataAction(completionEvents ++ failureEvents)
  }

  private def randomNumberString: String = Random.nextInt.toString.stripPrefix("-")

  private def completedJobMetadataEvents(jobKey: JobKey, executionStatus: ExecutionStatus, returnCode: Option[Int]) = {
    val returnCodeEvent = returnCode map { rc =>
      List(MetadataEvent(metadataKey(jobKey, CallMetadataKeys.ReturnCode), MetadataValue(rc)))
    }

    List(
      MetadataEvent(metadataKey(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(executionStatus)),
      MetadataEvent(metadataKey(jobKey, CallMetadataKeys.End), MetadataValue(OffsetDateTime.now))
    ) ++ returnCodeEvent.getOrElse(List.empty)
  }

  /**
    * Attempt to start all runnable jobs and return updated state data.  This will create a new copy
    * of the state data including new pending persists.
    */
  @tailrec
  private def startRunnableScopes(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {
    val runnableScopes = data.executionStore.runnableScopes.toList
    val runnableCalls = runnableScopes.view collect { case k if k.scope.isInstanceOf[Call] => k } sortBy { k =>
      (k.scope.fullyQualifiedName, k.index.getOrElse(-1)) } map { _.tag }
    if (runnableCalls.nonEmpty) workflowLogger.info(s"Starting calls: " + runnableCalls.mkString(", "))

    // Each process returns a Try[WorkflowExecutionDiff], which, upon success, contains potential changes to be made to the execution store.
    val executionDiffs = runnableScopes map {
      case k: BackendJobDescriptorKey => processRunnableJob(k, data)
      case k: ScatterKey => processRunnableScatter(k, data)
      case k: CollectorKey => processRunnableCollector(k, data)
      case k =>
        val exception = new UnsupportedOperationException(s"Unknown entry in execution store: ${k.tag}")
        self ! JobInitializationFailed(k, exception)
        Failure(exception)
    }

    TryUtil.sequence(executionDiffs) match {
      case Success(diffs) if diffs.exists(_.containsNewEntry) => startRunnableScopes(data.mergeExecutionDiffs(diffs))
      case Success(diffs) => data.mergeExecutionDiffs(diffs)
      case Failure(e) => data
    }
  }

  private def pushNewJobMetadata(jobKey: BackendJobDescriptorKey, backendName: String) = {
    val startEvents = List(
      MetadataEvent(metadataKey(jobKey, CallMetadataKeys.Start), MetadataValue(OffsetDateTime.now)),
      MetadataEvent(metadataKey(jobKey, CallMetadataKeys.Backend), MetadataValue(backendName))
    )

    serviceRegistryActor ! PutMetadataAction(startEvents)
  }

  private def pushPreparedJobMetadata(jobKey: BackendJobDescriptorKey, jobInputs: Map[LocallyQualifiedName, WdlValue]) = {
    val inputEvents = jobInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKey(jobKey, s"${CallMetadataKeys.Inputs}")))
      case inputs =>
        inputs flatMap {
          case (inputName, inputValue) =>
            wdlValueToMetadataEvents(metadataKey(jobKey, s"${CallMetadataKeys.Inputs}:$inputName"), inputValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(inputEvents)
  }

  private def processRunnableJob(jobKey: BackendJobDescriptorKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        val exception = new IllegalStateException(s"$tag $message")
        workflowLogger.error(exception, s"$tag $message")
        throw exception
      case Some(backendName) =>
        factories.get(backendName) match {
          case Some(factory) =>
            val jobPreparationActorName = s"${workflowDescriptor.id}-BackendPreparationActor-${jobKey.tag}"
            val jobPreparationActor = context.actorOf(JobPreparationActor.props(data, jobKey, factory), jobPreparationActorName)
            pushNewJobMetadata(jobKey, backendName)
            jobPreparationActor ! JobPreparationActor.Start
            Success(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Starting)))
          case None =>
            throw new WorkflowExecutionException(List(new Exception(s"Could not get BackendLifecycleActor for backend $backendName")))
        }
    }
  }

  private def processRunnableScatter(scatterKey: ScatterKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    val lookup = data.hierarchicalLookup(scatterKey.scope, None) _

    scatterKey.scope.collection.evaluate(lookup, data.expressionLanguageFunctions) map {
      case a: WdlArray => WorkflowExecutionDiff(scatterKey.populate(a.value.size) + (scatterKey -> ExecutionStatus.Done))
      case v: WdlValue => throw new Throwable("Scatter collection must evaluate to an array")
    }
  }

  private def processRunnableCollector(collector: CollectorKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    val shards = data.executionStore.findShardEntries(collector) collect { case (k: BackendJobDescriptorKey, v) if v == ExecutionStatus.Done => k }
    data.outputStore.generateCollectorOutput(collector, shards) match {
      case Failure(e) => Failure(new RuntimeException(s"Failed to collect output shards for call ${collector.tag}"))
      case Success(outputs) => self ! ScatterCollectionSucceededResponse(collector, outputs)
        Success(WorkflowExecutionDiff(Map(collector -> ExecutionStatus.Starting)))
    }
  }

  private def metadataKey(jobKey: JobKey, myKey: String) = MetadataKey(workflowDescriptor.id, Option(MetadataJobKey(jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)), myKey)
}
