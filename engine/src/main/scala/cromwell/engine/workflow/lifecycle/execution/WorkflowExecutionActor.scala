package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, JobKey}
import cromwell.core.{WorkflowId, _}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.NotStarted
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.WorkflowExecutionActorState
import cromwell.engine.{EngineWorkflowDescriptor, ExecutionStatus, workflow}
import cromwell.services.MetadataServiceActor._
import lenthall.exception.ThrowableAggregation
import wdl4s._
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlValue}

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowExecutionActor {

  /**
    * States
    */
  sealed trait WorkflowExecutionActorState { def terminal = false }
  sealed trait WorkflowExecutionActorTerminalState extends WorkflowExecutionActorState { override val terminal = true }

  case object WorkflowExecutionPendingState extends WorkflowExecutionActorState
  case object WorkflowExecutionInProgressState extends WorkflowExecutionActorState
  case object WorkflowExecutionSuccessfulState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionFailedState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionAbortedState extends WorkflowExecutionActorTerminalState

  /**
    * Commands
    */
  sealed trait WorkflowExecutionActorCommand
  case object StartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object RestartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object AbortExecutingWorkflowCommand extends WorkflowExecutionActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowExecutionActorResponse
  case object WorkflowExecutionSucceededResponse extends WorkflowExecutionActorResponse
  case object WorkflowExecutionAbortedResponse extends WorkflowExecutionActorResponse
  final case class WorkflowExecutionFailedResponse(reasons: Seq[Throwable]) extends WorkflowExecutionActorResponse

  /**
    * Internal control flow messages
    */
  private case class JobInitializationFailed(throwable: Throwable)
  private case class ScatterCollectionFailedResponse(collectorKey: CollectorKey, throwable: Throwable)
  private case class ScatterCollectionSucceededResponse(collectorKey: CollectorKey, outputs: CallOutputs)

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
}

final case class WorkflowExecutionActor(workflowId: WorkflowId,
                                        workflowDescriptor: EngineWorkflowDescriptor,
                                        serviceRegistryActor: ActorRef)
  extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] {

  import WorkflowExecutionActor._
  import lenthall.config.ScalaConfig._

  val tag = self.path.name
  private lazy val DefaultMaxRetriesFallbackValue = 10

  implicit val actorSystem = context.system
  implicit val ec = context.dispatcher

  // TODO: We should probably create a trait which loads all the configuration (once per application), and let classes mix it in
  // to avoid doing ConfigFactory.load() at multiple places
  val MaxRetries = ConfigFactory.load().getIntOption("system.max-retries") match {
    case Some(value) => value
    case None =>
      log.warning(s"Failed to load the max-retries value from the configuration. Defaulting back to a value of `$DefaultMaxRetriesFallbackValue`.")
      DefaultMaxRetriesFallbackValue
  }

  private val factories = TryUtil.sequenceMap(workflowDescriptor.backendAssignments.values.toSet[String] map { backendName =>
    backendName -> CromwellBackends.shadowBackendLifecycleFactory(backendName)
  } toMap) recover {
    case e => throw new RuntimeException("Could not instantiate backend factories", e)
  } get

  private val configs = TryUtil.sequenceMap(workflowDescriptor.backendAssignments.values.toSet[String] map { backendName =>
    backendName -> BackendConfiguration.backendConfigurationDescriptor(backendName)
  } toMap) recover {
    case e => throw new RuntimeException("Could not instantiate backend configurations", e)
  } get

  // Initialize the StateData with ExecutionStore (all calls as NotStarted) and SymbolStore
  startWith(
    WorkflowExecutionPendingState,
    WorkflowExecutionActorData(
      workflowDescriptor,
      executionStore = buildInitialExecutionStore(),
      outputStore = OutputStore.empty))

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
    case Event(AbortExecutingWorkflowCommand, _) =>
      context.parent ! WorkflowExecutionAbortedResponse
      goto(WorkflowExecutionAbortedState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, actorProps), stateData) =>
      context.actorOf(actorProps, buildJobExecutionActorName(jobDescriptor)) ! ExecuteJobCommand
      stay() using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobDescriptor.key -> ExecutionStatus.Running)))
    case Event(BackendJobPreparationFailed(jobKey, t), stateData) =>
      log.error(s"Failed to start job $jobKey", t)
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
    case Event(SucceededResponse(jobKey, callOutputs), stateData) =>
      handleCallSuccessful(jobKey, callOutputs, stateData)
    case Event(FailedNonRetryableResponse(jobKey, reason, _), stateData) =>
      log.warning(s"Job ${jobKey.call.fullyQualifiedName} failed! Reason: ${reason.getMessage}", reason)
      context.parent ! WorkflowExecutionFailedResponse(List(reason))
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
    case Event(FailedRetryableResponse(jobKey, reason, _), stateData) =>
      log.warning(s"Job ${jobKey.tag} failed with a retryable failure: ${reason.getMessage}")
      handleRetryableFailure(jobKey)
    case Event(JobInitializationFailed(reason), stateData) =>
      log.warning(s"Jobs failed to initialize: $reason")
      goto(WorkflowExecutionFailedState)
    case Event(ScatterCollectionSucceededResponse(jobKey, callOutputs), stateData) =>
      handleCallSuccessful(jobKey, callOutputs, stateData)
    case Event(AbortExecutingWorkflowCommand, stateData) => ??? // TODO: Implement!
    case Event(_, _) => ??? // TODO: Lots of extra stuff to include here...
  }

  when(WorkflowExecutionSuccessfulState) {
    FSM.NullFunction
  }
  when(WorkflowExecutionFailedState) {
    FSM.NullFunction
  }
  when(WorkflowExecutionAbortedState) {
    FSM.NullFunction
  }

  whenUnhandled {
    case Event(MetadataPutFailed(action, error), _) =>
      // Do something useful here??
      log.warning(s"$tag Put failed for Metadata action $action : ${error.getMessage}")
      stay
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage in state: $stateName")
      stay
  }

  onTransition {
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState.")
  }

  onTransition {
    case _ -> WorkflowExecutionSuccessfulState =>
      pushOutputsToMetadataService(nextStateData)
      context.parent ! WorkflowExecutionSucceededResponse
    case _ -> WorkflowExecutionAbortedState => context.parent ! WorkflowExecutionAbortedResponse
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.info(s"$tag done. Shutting down.")
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey) = {
    // We start with index 1 for #attempts, hence invariant breaks only if jobKey.attempt > MaxRetries
    if (jobKey.attempt <= MaxRetries) {
      val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
      log.info(s"Retrying job execution for ${newJobKey.tag}")
      /** Currently, we update the status of the old key to Preempted, and add a new entry (with the #attempts incremented by 1)
        * to the execution store with status as NotStarted. This allows startRunnableCalls to re-execute this job */
      val executionDiff = WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Preempted, newJobKey -> ExecutionStatus.NotStarted))
      val newData = stateData.mergeExecutionDiff(executionDiff)
      stay() using startRunnableScopes(newData)
    } else {
      log.warning(s"Exhausted maximum number of retries for job ${jobKey.tag}. Failing.")
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))
    }
  }

  def handleCallSuccessful(jobKey: JobKey, outputs: CallOutputs, data: WorkflowExecutionActorData) = {
    log.info(s"Job ${jobKey.tag} succeeded! Outputs: ${outputs.mkString("\n")}")
    val newData = data.jobExecutionSuccess(jobKey, outputs)

    if (newData.isWorkflowComplete) {
      log.info(newData.outputsJson())
      goto(WorkflowExecutionSuccessfulState) using newData
    }
    else
      stay() using startRunnableScopes(newData)
  }

  /**
    * Attempt to start all runnable jobs and return updated state data.  This will create a new copy
    * of the state data including new pending persists.
    */
  @tailrec
  private def startRunnableScopes(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {
    val runnableScopes = data.executionStore.runnableScopes.toList
    val runnableCalls = runnableScopes.view collect { case k if k.scope.isInstanceOf[Call] => k } sortBy { _.index.getOrElse(-1) } map { _.tag }
    if (runnableCalls.nonEmpty) log.info(s"Starting calls: " + runnableCalls.mkString(", "))

    // Each process*** returns a Try[WorkflowExecutionDiff], which, upon success, contains potential changes to be made to the execution store.
    val executionDiffs = runnableScopes map {
      case k: BackendJobDescriptorKey => processRunnableJob(k, data)
      case k: ScatterKey => processRunnableScatter(k, data)
      case k: CollectorKey => processRunnableCollector(k, data)
      case k => Failure(new UnsupportedOperationException(s"Unknown entry in execution store: ${k.tag}"))
    }

    TryUtil.sequence(executionDiffs) match {
      case Success(diffs) if diffs.exists(_.containsNewEntry) => startRunnableScopes(data.mergeExecutionDiffs(diffs))
      case Success(diffs) => data.mergeExecutionDiffs(diffs)
      case Failure(e) =>
        self ! JobInitializationFailed(e)
        data
    }
  }

  private def processRunnableJob(jobKey: BackendJobDescriptorKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        log.error(s"$tag $message")
        throw new IllegalStateException(s"$tag $message")
      case Some(backendName) =>
        (configs.get(backendName), factories.get(backendName)) match {
          case (Some(configDescriptor), Some(factory)) =>
            val jobPreparationActorName = s"${workflowDescriptor.id}-BackendPreparationActor-${jobKey.tag}"
            val jobPreparationActor = context.actorOf(JobPreparationActor.props(data, jobKey, factory, configDescriptor), jobPreparationActorName)
            jobPreparationActor ! JobPreparationActor.Start
            Success(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Starting)))
          case (c, f) =>
            val noConf = if (c.isDefined) None else Option(new Exception(s"Could not get BackendConfigurationDescriptor for backend $backendName"))
            val noFactory = if (f.isDefined) None else Option(new Exception(s"Could not get BackendLifecycleActor for backend $backendName"))
            val errors = List(noConf, noFactory).flatten
            errors foreach(error => log.error(error.getMessage, error))
            throw new WorkflowExecutionException(errors)
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

  private def pushOutputsToMetadataService(data: WorkflowExecutionActorData): Unit = {
    import workflow._
    val keyValues = data.outputStore.store.flatMap {
      case (key, value) => value map (entry => s"${key.call.fullyQualifiedName}.${entry.name}" -> entry.wdlValue.map(_.toWdlString).getOrElse("NA"))
    }
    val metadataEventMsgs = keyValues map{ case (k,v) => MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Outputs}:$k"), MetadataValue(v))}
    metadataEventMsgs foreach ( serviceRegistryActor ! PutMetadataAction(_) )
  }

}
