package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, BackendJobExecutionFailedRetryableResponse, BackendJobExecutionSucceededResponse, ExecuteJobCommand}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, JobKey}
import cromwell.core.{WorkflowId, _}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.{ExecutionStatus, NotStarted}
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.JobInputEvaluator
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.{WorkflowExecutionActorData, WorkflowExecutionActorState}
import cromwell.engine.{EngineWorkflowDescriptor, ExecutionStatus}
import cromwell.webservice.WdlValueJsonFormatter
import lenthall.exception.ThrowableAggregation
import wdl4s._
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

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

  type ExecutionStore = Map[JobKey, ExecutionStatus]
  type ExecutionStoreEntry = (JobKey, ExecutionStatus)

  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputCallKey(call: Call, index: ExecutionIndex)
  type OutputStore = Map[OutputCallKey, Traversable[OutputEntry]]

  /**
    * State data
    */
  final case class WorkflowExecutionActorData(executionStore: ExecutionStore,
                                              outputStore: OutputStore) {

    /** This method updates: The ExecutionStore with the updated status and the symbol cache with the new outputs */
    def updateJob(jobKey: BackendJobDescriptorKey,
                  outputs: CallOutputs,
                  status: ExecutionStatus): WorkflowExecutionActorData = {
      this.copy(executionStore = executionStore + (jobKey -> status),
        outputStore = outputStore ++ updateSymbolStoreEntry(jobKey, outputs))
    }

    /** Add the outputs for the specified `JobKey` to the symbol cache. */
    private def updateSymbolStoreEntry(jobKey: BackendJobDescriptorKey, outputs: CallOutputs) = {
      val newOutputEntries = outputs map {
        case (name, value) => OutputEntry(name, value.wdlValue.wdlType, Option(value.wdlValue))
      }

      outputStore + (OutputCallKey(jobKey.call, jobKey.index) -> newOutputEntries)
    }

    /** Checks if the workflow is completed by scanning through the executionStore */
    def isWorkflowComplete: Boolean = {
      def isDone(executionStatus: ExecutionStatus): Boolean = (executionStatus == ExecutionStatus.Done) || (executionStatus == ExecutionStatus.Preempted)
      executionStore.values.forall(isDone)
    }

    def containsFailedJob: Boolean = {
      executionStore.values.exists(_ == ExecutionStatus.Failed)
    }

    /** Updates the status of a a job by the new entry */
    def updateExecutionStoreStatus(entry: (JobKey, ExecutionStatus)): WorkflowExecutionActorData = this.copy(executionStore = executionStore + entry)
  }

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

  private case class JobInitializationFailed(jobKey: JobKey, throwable: Throwable)

  case class WorkflowExecutionException(override val throwables: List[Throwable]) extends ThrowableAggregation {
    override val exceptionContext = s"WorkflowExecutionActor"
  }

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(WorkflowExecutionActor(workflowId, workflowDescriptor))
}

final case class WorkflowExecutionActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] {

  import WorkflowExecutionActor._
  import lenthall.config.ScalaConfig._

  val tag = self.path.name
  private lazy val DefaultMaxRetriesFallbackValue = 10

  private val calls = workflowDescriptor.backendDescriptor.workflowNamespace.workflow.calls
  private val inputResolver = new JobInputEvaluator(workflowDescriptor)

  // TODO: We should probably create a trait which loads all the configuration (once per application), and let classes mix it in
  // to avoid doing ConfigFactory.load() at multiple places
  val MaxRetries = ConfigFactory.load().getIntOption("system.max-retries") match {
    case Some(value) => value
    case None =>
      log.warning(s"Failed to load the max-retries value from the configuration. Defaulting back to a value of `$DefaultMaxRetriesFallbackValue`.")
      DefaultMaxRetriesFallbackValue
  }

  // Initialize the StateData with ExecutionStore (all calls as NotStarted) and SymbolStore
  startWith(
    WorkflowExecutionPendingState,
    WorkflowExecutionActorData(
      executionStore = (calls map (BackendJobDescriptorKey(_, None, 1) -> NotStarted)) toMap,
      outputStore = Map.empty))

  // Return a more meaningful value that execution status
  private def executeJob(jobKey: BackendJobDescriptorKey,
                         outputStore: OutputStore,
                         configDescriptor: BackendConfigurationDescriptor,
                         factory: BackendLifecycleActorFactory): Try[ExecutionStatus] = {
    // FIXME This is a potential bottleneck as input evaluation can be arbitrarily long and executes code coming from the backend.
    // It should probably be actorified / isolated
    inputResolver.resolveAndEvaluate(jobKey, factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, jobKey, configDescriptor), outputStore) map { inputs =>
      val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, inputs)
      val actorName = s"${jobDescriptor.descriptor.id}-BackendExecutionActor-${jobKey.tag}"
      val jobExecutionActor = context.actorOf(
        factory.jobExecutionActorProps(
          jobDescriptor,
          BackendConfigurationDescriptor(configDescriptor.backendConfig, configDescriptor.globalConfig)
        ),
        actorName
      )
      jobExecutionActor ! ExecuteJobCommand
      ExecutionStatus.Starting
    }
  }

  private def startJob(jobKey: BackendJobDescriptorKey, symbolCache: OutputStore): Try[ExecutionStatus] = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        log.error(s"$tag $message")
        throw new IllegalStateException(s"$tag $message")
      case Some(backendName) =>
        // TODO these shouldn't be re-instantiated for every call
        val attemptedConfigurationDescriptor = BackendConfiguration.backendConfigurationDescriptor(backendName)
        val attemptedActorFactory = CromwellBackends.shadowBackendLifecycleFactory(backendName)

        (attemptedConfigurationDescriptor, attemptedActorFactory) match {
          case (Success(configDescriptor), Success(factory)) =>
            executeJob(jobKey, symbolCache, configDescriptor, factory)
          case (_, _) =>
            val errors = List(
              attemptedActorFactory.failed.map(new Exception(s"Could not get BackendLifecycleActor for backend $backendName", _)).toOption,
              attemptedConfigurationDescriptor.failed.map(new Exception(s"Could not get BackendConfigurationDescriptor for backend $backendName", _)).toOption
            ).flatten
            errors foreach(error => log.error(error.getMessage, error))
            throw new WorkflowExecutionException(errors)
        }
    }
  }

  when(WorkflowExecutionPendingState) {
    case Event(StartExecutingWorkflowCommand, stateData) =>
      val data = startRunnableJobs(stateData)
      goto(WorkflowExecutionInProgressState) using data
    case Event(RestartExecutingWorkflowCommand, _) =>
      // TODO: Restart executing
      goto(WorkflowExecutionInProgressState)
    case Event(AbortExecutingWorkflowCommand, _) =>
      context.parent ! WorkflowExecutionAbortedResponse
      goto(WorkflowExecutionAbortedState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(BackendJobExecutionSucceededResponse(jobKey, callOutputs), stateData) =>
      log.info(s"Job ${jobKey.call.fullyQualifiedName} succeeded! Outputs: ${callOutputs.mkString("\n")}")
      val newData = stateData.updateJob(jobKey, callOutputs, ExecutionStatus.Done)
      if (newData.isWorkflowComplete) {
        printOutputs(newData.outputStore)
        goto(WorkflowExecutionSuccessfulState) using newData
      }
      else
        stay() using startRunnableJobs(newData)
    case Event(BackendJobExecutionFailedResponse(jobKey, reason), stateData) =>
      log.warning(s"Job ${jobKey.call.fullyQualifiedName} failed! Reason: ${reason.getMessage}", reason)
      goto(WorkflowExecutionFailedState) using stateData.updateExecutionStoreStatus(jobKey -> ExecutionStatus.Failed)
    case Event(BackendJobExecutionFailedRetryableResponse(jobKey, reason), stateData) =>
      log.warning(s"Job ${jobKey.tag} failed with a retryable failure: ${reason.getMessage}")
      handleRetryableFailure(jobKey)
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      log.warning(s"Job ${jobKey.tag} failed to initialize: $reason")
      goto(WorkflowExecutionFailedState)
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
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage in state: $stateName")
      stay
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.info(s"$tag done. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState.")
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey) = {
    // We start with index 1 for #attempts, hence invariant breaks only if jobKey.attempt > MaxRetries
    if (jobKey.attempt <= MaxRetries) {
      val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
      log.info(s"Retrying job execution for ${newJobKey.tag}")
      /** Currently, we update the status of the old key to Preempted, and add a new entry (with the #attempts incremented by 1)
        * to the execution store with status as NotStarted. This allows startRunnableCalls to re-execute this job */
      val newData = stateData.updateExecutionStoreStatus(jobKey, ExecutionStatus.Preempted).updateExecutionStoreStatus(newJobKey -> ExecutionStatus.NotStarted)
      stay() using startRunnableJobs(newData)
    } else {
      log.warning(s"Exhausted maximum number of retries for job ${jobKey.tag}. Failing.")
      goto(WorkflowExecutionFailedState) using stateData.updateExecutionStoreStatus(jobKey -> ExecutionStatus.Failed)
    }
  }

  private def upstreamEntries(entry: JobKey, prerequisiteScope: Scope, executionStore: ExecutionStore): Seq[ExecutionStoreEntry] = {
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
        executionStore filter { case (k, _) => k.scope == prerequisiteScope && k.index == entry.index } toSeq

      /**
        * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
        * on every shard of the pre-requisite scope to finish.
        */
      case _ =>
        executionStore filter { case (k, _) => k.scope == prerequisiteScope && k.index.isEmpty } toSeq
    }
  }

  private def arePrerequisitesDone(key: JobKey, executionStore: ExecutionStore): Boolean = {
    def isDone(e: JobKey): Boolean =
      executionStore exists { case (k, s) => k.scope == e.scope && k.index == e.index && s == ExecutionStatus.Done }

    val upstream = key.scope.prerequisiteScopes.map(s => upstreamEntries(key, s, executionStore))
    // TODO: Check downstream for Scatter calls (?)
    executionStore.filter(upstream.flatten) forall (x => isDone(x._1))
  }

  private def isRunnable(entry: ExecutionStoreEntry, executionStore: ExecutionStore) = {
    entry match {
      case (key, NotStarted) => arePrerequisitesDone(key, executionStore)
      case _ => false
    }
  }

  /**
    * Attempt to start all runnable jobs and return updated state data.  This will create a new copy
    * of the state data including new pending persists.
    */
  @tailrec
  private def startRunnableJobs(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {

    val executionStore = data.executionStore

    val runnableEntries = executionStore filter(isRunnable(_, executionStore))
    val runnableCalls = runnableEntries collect { case (k: BackendJobDescriptorKey, v) => k.scope }
    if (runnableCalls.nonEmpty)
      log.info(s"Starting calls: " + runnableCalls.map(_.fullyQualifiedName).toSeq.sorted.mkString(", "))

    val entries: Map[JobKey, Try[ExecutionStatus]] = runnableEntries map {
      case (k: BackendJobDescriptorKey, _) => k -> startJob(k, data.outputStore)
      case (k, v) =>
        val message = s"Unknown entry in execution store:\nKEY: ${k.tag}\nVALUE:$v"
        log.error(message)
        k -> Failure(new UnsupportedOperationException(message))
    }

    val updatedData = entries.foldLeft(data)((newData, entry) => {
      entry match {
        case (key, status) if status.isSuccess => newData.updateExecutionStoreStatus(key, ExecutionStatus.Starting)
        case (key, status) =>
          self ! JobInitializationFailed(key, status.failed.get)
          newData.updateExecutionStoreStatus(key, ExecutionStatus.Failed)
      }
    })

    if (entries.nonEmpty) startRunnableJobs(updatedData) else updatedData
  }

  private def printOutputs(outputStore: OutputStore) = {
    // Printing the final outputs, temporarily here until SingleWorkflowManagerActor is made in-sync with the shadow mode
    import WdlValueJsonFormatter._
    import spray.json._
    val workflowOutputs = outputStore flatMap {
      case (key, outputs) => outputs map { output =>
        s"${key.call.fullyQualifiedName}.${output.name}" -> (output.wdlValue map { _.valueString } getOrElse "N/A")
      }
    }
    log.info("Workflow complete. Final Outputs: \n" + workflowOutputs.toJson.prettyPrint)
  }
}
