package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, BackendJobExecutionSucceededResponse, ExecuteJobCommand}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, JobKey}
import cromwell.core.{WorkflowId, _}
import cromwell.engine.ExecutionStatus._
import cromwell.engine.{EngineWorkflowDescriptor, ExecutionStatus, _}
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor._
import cromwell.webservice.WdlValueJsonFormatter
import lenthall.exception.ThrowableAggregation
import wdl4s._
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

  case class SymbolCacheKey(scopeName: String, input: Boolean)

  type SymbolCache = Map[SymbolCacheKey, Traversable[SymbolStoreEntry]]
  type ExecutionStore = Map[JobKey, ExecutionStatus]
  type ExecutionStoreEntry = (JobKey, ExecutionStatus)

  /**
    * State data
    */
  final case class WorkflowExecutionActorData(executionStore: ExecutionStore,
                                              symbolCache: SymbolCache) {

    /** This method updates: The ExecutionStore with the updated status, the symbol cache with the new outputs, and appends to the overall CallOutputs */
    def updateJob(jobKey: JobKey,
                  outputs: CallOutputs,
                  status: ExecutionStatus): WorkflowExecutionActorData = {
      this.copy(executionStore = executionStore + (jobKey -> status),
        symbolCache = symbolCache ++ updateSymbolStoreEntry(jobKey, outputs))
    }

    /** Add the outputs for the specified `JobKey` to the symbol cache. */
    private def updateSymbolStoreEntry(jobKey: JobKey, outputs: CallOutputs) = {
      val newEntriesMap = outputs map { case (lqn, value) =>
        val storeKey = SymbolStoreKey(jobKey.scope.fullyQualifiedName, lqn, jobKey.index, input = false)
        // TODO: SymbolStoreEntry should no longer contain the symbol hashes
        new SymbolStoreEntry(storeKey, value.wdlValue.wdlType, Option(value.wdlValue), None)
      } groupBy { entry => SymbolCacheKey(entry.scope, entry.isInput) }

      newEntriesMap map { case (key, entries) =>
        // SymbolCache is essentially a MultiMap, but that's a trait only for mutable Maps.
        key -> (symbolCache.getOrElse(key, Seq.empty) ++ entries)
      }
    }

    /** Checks if the workflow is completed by scanning through the executionStore */
    def isWorkflowComplete: Boolean = {
      def isDone(executionStatus: ExecutionStatus): Boolean = executionStatus == ExecutionStatus.Done
      executionStore.values.forall(isDone)
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

  case class WorkflowExecutionException(override val throwables: List[Throwable]) extends ThrowableAggregation {
    override val exceptionContext = s"WorkflowExecutionActor"
  }

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(WorkflowExecutionActor(workflowId, workflowDescriptor))
}

final case class WorkflowExecutionActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] {

  import WorkflowExecutionActor._

  val tag = self.path.name

  private val calls = workflowDescriptor.backendDescriptor.workflowNamespace.workflow.calls

  // Initialize the StateData with ExecutionStore (all calls as NotStarted) and SymbolStore
  startWith(
    WorkflowExecutionPendingState,
    WorkflowExecutionActorData(
      executionStore = (calls map (BackendJobDescriptorKey(_, None, 1) -> NotStarted)) toMap,
      symbolCache = buildSymbolStoreEntries.groupBy(entry => SymbolCacheKey(entry.scope, entry.isInput))))

  /** PBE: the return value of WorkflowExecutionActorState is just temporary.
    *      This should probably return a Try[BackendJobDescriptor], Unit, Boolean,
    *      Try[ActorRef], or something to indicate if the job was started
    *      successfully.  Or, if it can fail to start, some indication of why it
    *      failed to start
    */
  private def startJob(jobKey: BackendJobDescriptorKey,
                       symbolCache: SymbolCache,
                       configDescriptor: BackendConfigurationDescriptor,
                       factory: BackendLifecycleActorFactory): Try[ExecutionStatus] = {

      val call = jobKey.call
      fetchLocallyQualifiedInputs(jobKey, symbolCache) match {
        case Success(symbolStoreForCall) =>
          val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, symbolStoreForCall)
          val actorName = s"${jobDescriptor.descriptor.id}-BackendExecutionActor-${jobKey.call.fullyQualifiedName}-${jobKey.index.getOrElse("")}-${jobKey.attempt}"
          val jobExecutionActor = context.actorOf(
            factory.jobExecutionActorProps(
              jobDescriptor,
              BackendConfigurationDescriptor(configDescriptor.backendConfig, configDescriptor.globalConfig)
            ),
            actorName
          )
          jobExecutionActor ! ExecuteJobCommand
          Success(ExecutionStatus.Starting)
        case Failure(reason) =>
          log.error(s"Failed to fetch locally qualified inputs for call ${call.fullyQualifiedName}", reason)
          throw new WorkflowExecutionException(List(reason))
      }
  }

  private def startJob(jobKey: BackendJobDescriptorKey, symbolCache: SymbolCache): Try[ExecutionStatus] = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        log.error(s"$tag $message")
        throw new IllegalStateException(s"$tag $message")
      case Some(backendName) =>
        val attemptedConfigurationDescriptor = BackendConfiguration.backendConfigurationDescriptor(backendName)
        val attemptedActorFactory = CromwellBackends.shadowBackendLifecycleFactory(backendName)

        (attemptedConfigurationDescriptor, attemptedActorFactory) match {
          case (Success(configDescriptor), Success(factory)) =>
            startJob(jobKey, symbolCache, configDescriptor, factory)
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
        printOutputs(newData)
        goto(WorkflowExecutionSuccessfulState) using newData
      }
      else
        stay() using startRunnableJobs(newData)
    case Event(BackendJobExecutionFailedResponse(jobKey, reason), stateData) =>
      log.warning(s"Job ${jobKey.call.fullyQualifiedName} failed! Reason: $reason")
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

  private def fetchCallInputEntries(callKey: JobKey, symbolCache: SymbolCache): Traversable[SymbolStoreEntry] = {
    symbolCache.getOrElse(SymbolCacheKey(callKey.scope.fullyQualifiedName, input = true), Seq.empty)
  }

  private def buildSymbolStoreEntries: Traversable[SymbolStoreEntry] = {
    val actualInputs = workflowDescriptor.backendDescriptor.inputs ++ workflowDescriptor.declarations
    val inputSymbols = actualInputs map {
      case (name, value) => SymbolStoreEntry(name, value, None, input = true)
    }

    val callSymbols = for {
      call <- workflowDescriptor.namespace.workflow.calls
      (k, v) <- call.inputMappings
    } yield SymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, None, input = true)

    inputSymbols.toSet ++ callSymbols.toSet
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
      case (key, ExecutionStatus.NotStarted) => arePrerequisitesDone(key, executionStore)
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
      case (k: BackendJobDescriptorKey, _) => k -> startJob(k, data.symbolCache)
      case (k, v) =>
        val message = s"Unknown entry in execution store:\nKEY: ${k.tag}\nVALUE:$v"
        log.error(message)
        k -> Failure(new UnsupportedOperationException(message))
    }

    entries.filter(_._2.isFailure) foreach { case (key, status) => log.error(s"Failed to start Job ${key.scope.fullyQualifiedName}: ${status.failed.get.getMessage}")}

    val updatedData = entries.filter(_._2.isSuccess).foldLeft(data)((newData, entry) => newData.updateExecutionStoreStatus(entry._1, entry._2.get))
    if (entries.nonEmpty) startRunnableJobs(updatedData) else updatedData
  }

  private def fetchLocallyQualifiedInputs(callKey: JobKey, symbolCache: SymbolCache): Try[Map[String, WdlValue]] = Try {
    val entries = fetchCallInputEntries(callKey, symbolCache)
    entries.map { entry =>
      val value = entry.wdlValue match {
        case Some(v) => v
        case _ => throw new WdlExpressionException("Unknown error")
      }

      // TODO: Coercion to happen here? We don't have EngineFunctions here because of which the abpve pattern match
      // cannot have the WdlExpressions evaluated.
      // val coercedValue = value.flatMap(x => declaration.wdlType.coerceRawValue(x))
      entry.key.name -> value
    }.toMap
  }

  private def printOutputs(stateData: WorkflowExecutionActorData) = {
    // Printing the final outputs, temporarily here until SingleWorkflowManagerActor is made in-sync with the shadow mode
    import WdlValueJsonFormatter._
    import spray.json._
    val workflowOutputs = stateData.symbolCache.flatMap(_._2).collect {
      case x if x.isOutput => s"${x.key.scope}.${x.key.name}" -> x.wdlValue
    }.toMap
    log.info("Workflow complete. Final Outputs: \n" + workflowOutputs.toJson.prettyPrint)
  }
}
