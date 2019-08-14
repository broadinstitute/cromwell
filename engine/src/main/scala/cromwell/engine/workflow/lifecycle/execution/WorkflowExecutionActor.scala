package cromwell.engine.workflow.lifecycle.execution

import java.util.concurrent.atomic.AtomicInteger

import _root_.wdl.draft2.model._
import akka.actor.{Scope => _, _}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.Checked
import common.exception.{AggregatedException, AggregatedMessageException, MessageAggregation}
import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked._
import common.validation.Validation._
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend._
import cromwell.backend.standard.callcaching.BlacklistCache
import cromwell.core.Dispatcher._
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.core.io.AsyncIo
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData.DataStoreUpdate
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor
import cromwell.engine.workflow.lifecycle.execution.keys.ExpressionKey.{ExpressionEvaluationFailedResponse, ExpressionEvaluationSucceededResponse}
import cromwell.engine.workflow.lifecycle.execution.keys._
import cromwell.engine.workflow.lifecycle.execution.stores.{ActiveExecutionStore, ExecutionStore}
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}
import cromwell.engine.workflow.workflowstore.{RestartableAborting, StartableState}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{CallMetadataKeys, MetadataEvent, MetadataValue}
import cromwell.util.StopAndLogSupervisor
import cromwell.webservice.EngineStatsActor
import net.ceedubs.ficus.Ficus._
import org.apache.commons.lang3.StringUtils
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.{ExposedExpressionNode, TaskCallInputExpressionNode}
import wom.values._

import scala.concurrent.duration._
import scala.language.postfixOps

case class WorkflowExecutionActor(params: WorkflowExecutionActorParams)
  extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] with WorkflowLogging with CallMetadataHelper with StopAndLogSupervisor with Timers {

  implicit val ec = context.dispatcher
  override val serviceRegistryActor = params.serviceRegistryActor
  val workflowDescriptor = params.workflowDescriptor
  override val workflowIdForLogging = workflowDescriptor.possiblyNotRootWorkflowId
  override val rootWorkflowIdForLogging = workflowDescriptor.rootWorkflowId
  override val workflowIdForCallMetadata = workflowDescriptor.id
  private val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)
  private val restarting = params.startState.restarted
  private val tag = s"WorkflowExecutionActor [UUID(${workflowDescriptor.id.shortString})]"

  private val DefaultTotalMaxJobsPerRootWf = 1000000
  private val DefaultMaxScatterSize = 1000000
  private val TotalMaxJobsPerRootWf = params.rootConfig.getOrElse("system.total-max-jobs-per-root-workflow", DefaultTotalMaxJobsPerRootWf)
  private val MaxScatterWidth = params.rootConfig.getOrElse("system.max-scatter-width-per-scatter", DefaultMaxScatterSize)

  private val backendFactories: Map[String, BackendLifecycleActorFactory] = {
    val factoriesValidation = workflowDescriptor.backendAssignments.values.toList
      .traverse[ErrorOr, (String, BackendLifecycleActorFactory)] {
      backendName => CromwellBackends.backendLifecycleFactoryActorByName(backendName) map { backendName -> _ }
    }

    factoriesValidation
      .map(_.toMap)
      .valueOr(errors => throw AggregatedMessageException("Could not instantiate backend factories", errors.toList))
  }


  val executionStore: ErrorOr[ActiveExecutionStore] = ExecutionStore(workflowDescriptor.callable, params.totalJobsByRootWf, TotalMaxJobsPerRootWf)

  // If executionStore returns a Failure about root workflow creating jobs more than total jobs per root workflow limit,
  // the WEA will fail by sending WorkflowExecutionFailedResponse to its parent and kill itself
  executionStore match {
    case Valid(validExecutionStore) =>
      startWith(
        WorkflowExecutionPendingState,
        WorkflowExecutionActorData(workflowDescriptor, ioEc, new AsyncIo(params.ioActor, GcsBatchCommandBuilder), params.totalJobsByRootWf, validExecutionStore)
      )
    case Invalid(e) =>
      val errorMsg = s"Failed to initialize WorkflowExecutionActor. Error: $e"
      workflowLogger.error(errorMsg)
      context.parent ! WorkflowExecutionFailedResponse(Map.empty, new Exception(errorMsg))

      // We start the actor with some state because if the actor is not started with a state before killing itself,
      // it throws NullPointerException as FSM.goto can't find the currentState
      startWith(
        WorkflowExecutionFailedState,
        WorkflowExecutionActorData(workflowDescriptor, ioEc, new AsyncIo(params.ioActor, GcsBatchCommandBuilder), params.totalJobsByRootWf, ExecutionStore.empty)
      )

      workflowLogger.debug("Actor failed to initialize. Stopping self.")
      context.stop(self)
  }

  private def sendHeartBeat(): Unit = timers.startSingleTimer(ExecutionHeartBeatKey, ExecutionHeartBeat, ExecutionHeartBeatInterval)

  when(WorkflowExecutionPendingState) {
    case Event(ExecuteWorkflowCommand, _) =>
      // Start HeartBeat
      sendHeartBeat()

      /*
        * Note that we don't record the fact that a workflow was failing, therefore we either restart in running or aborting state.
        * If the workflow was failing, it means at least one job failed in which case it'll still be failed when we get to it.
        * When that happens, we'll go to failing state.
        *
        * An effect of that is that up until the moment when we come across the failed job,
        * all backend jobs will be restarted with a Recover command which could potentially re-execute the job if the backend doesn't support
        * job recovery. A better way would be to record that the workflow was failing an restart in failing mode. However there will always be a gap
        * between when a job has failed and Cromwell is aware of it, or has time to persist that information.
       */
      params.startState match {
        case RestartableAborting => goto(WorkflowExecutionAbortingState)
        case _ => goto(WorkflowExecutionInProgressState)
      }
    case Event(EngineLifecycleActorAbortCommand, _) =>
      context.parent ! WorkflowExecutionAbortedResponse(Map.empty)
      goto(WorkflowExecutionAbortedState)
  }

  /* ********************* */
  /* ****** Running ****** */
  /* ********************* */

  when(WorkflowExecutionInProgressState) {
    // If we're done, the workflow is successful
    case Event(ExecutionHeartBeat, data) if data.done => handleWorkflowSuccessful(data)
  }

  /* ********************* */
  /* ****** Failing ****** */
  /* ********************* */

  when(WorkflowExecutionFailingState) {
    // If we're done, aggregate all job failures and fail the workflow
    case Event(ExecutionHeartBeat, data) if data.done =>
      val workflowFailure = AggregatedException("Workflow failed", data.jobFailures.values)
      context.parent ! WorkflowExecutionFailedResponse(data.jobExecutionMap, workflowFailure)
      goto(WorkflowExecutionFailedState)

    // A job not found here means we were trying to reconnect to a job that was likely never started. Indicate this in the message.
    case Event(JobFailedNonRetryableResponse(jobKey, _: JobNotFoundException, _), _) if restarting =>
      val benignException = new Exception("Cromwell server was restarted while this workflow was running. As part of the restart process, Cromwell attempted to reconnect to this job, however it was never started in the first place. This is a benign failure and not the cause of failure for this workflow, it can be safely ignored.")
      handleNonRetryableFailure(stateData, jobKey, benignException)
  }

  /* ********************** */
  /* ****** Aborting ****** */
  /* ********************** */

  when(WorkflowExecutionAbortingState) {
    // If we're done, the workflow is aborted
    case Event(ExecutionHeartBeat, data) if data.done =>
      context.parent ! WorkflowExecutionAbortedResponse(data.jobExecutionMap)
      goto(WorkflowExecutionAbortedState)

    case Event(JobAbortedResponse(jobKey), stateData) =>
      handleCallAborted(stateData, jobKey)

    case Event(SubWorkflowAbortedResponse(jobKey, executedKeys), stateData) =>
      handleCallAborted(stateData, jobKey, executedKeys)

    // Here we can't really know what the status of the job is. For now declare it aborted anyway but add some info in the metadata
    case Event(JobFailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, failure: JobReconnectionNotSupportedException, _), _) if restarting =>
      pushBackendStatusUnknown(jobKey)
      handleNonRetryableFailure(stateData, jobKey, failure, None)
  }

  when(WorkflowExecutionSuccessfulState) { FSM.NullFunction }
  when(WorkflowExecutionFailedState) { FSM.NullFunction }
  when(WorkflowExecutionAbortedState) { FSM.NullFunction }

  // Most of the Event handling is common to all states, so put it here. Specific behavior is added / overridden in each state
  whenUnhandled {
    case Event(ExecutionHeartBeat, data) if data.executionStore.needsUpdate =>
      val newData = startRunnableNodes(data)
      sendHeartBeat()
      stay() using newData
    case Event(ExecutionHeartBeat, data) if data.stalled =>
      handleWorkflowStalled(data)
    case Event(ExecutionHeartBeat, _) =>
      sendHeartBeat()
      stay()
    case Event(JobStarting(jobKey), stateData) =>
      pushStartingCallMetadata(jobKey)
      stay() using stateData
        .mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Starting)))
    case Event(JobRunning(key, inputs), stateData) =>
      pushRunningCallMetadata(key, inputs)
      stay() using stateData
        .mergeExecutionDiff(WorkflowExecutionDiff(Map(key -> ExecutionStatus.Running)))

    //Success
    // Job
    case Event(r: JobSucceededResponse, stateData) =>
      if (r.resultGenerationMode != RunOnBackend) {
        workflowLogger.info(s"Job results retrieved (${r.resultGenerationMode}): '${r.jobKey.call.fullyQualifiedName}' (scatter index: ${r.jobKey.index}, attempt ${r.jobKey.attempt})")
      }
      handleCallSuccessful(r.jobKey, r.jobOutputs, r.returnCode, stateData, Map.empty)
    // Sub Workflow
    case Event(SubWorkflowSucceededResponse(jobKey, descendantJobKeys, callOutputs), currentStateData) =>
      // Update call outputs to come from sub-workflow output ports:
      val subworkflowOutputs: Map[OutputPort, WomValue] = callOutputs.outputs flatMap { case (port, value) =>
        jobKey.node.subworkflowCallOutputPorts collectFirst {
          case p if p.outputToExpose.graphOutputPort eq port => p -> value
        }
      }

      if(subworkflowOutputs.size == callOutputs.outputs.size) {
        handleCallSuccessful(jobKey, CallOutputs(subworkflowOutputs), None, currentStateData, descendantJobKeys)
      } else {
        handleNonRetryableFailure(currentStateData, jobKey, new Exception(s"Subworkflow produced outputs: [${callOutputs.outputs.keys.mkString(", ")}], but we expected all of [${jobKey.node.subworkflowCallOutputPorts.map(_.internalName)}]"))
      }
    // Expression
    case Event(ExpressionEvaluationSucceededResponse(expressionKey, callOutputs), stateData) =>
      expressionKey.node match {
        case _: ExposedExpressionNode | _: ExpressionBasedGraphOutputNode =>
          workflowLogger.debug(s"Expression evaluation succeeded: '${expressionKey.node.fullyQualifiedName}' (scatter index: ${expressionKey.index}, attempt: ${expressionKey.attempt})")
        case _ => // No logging; anonymous node
      }

      handleDeclarationEvaluationSuccessful(expressionKey, callOutputs, stateData)

    // Failure
    // Initialization
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      handleNonRetryableFailure(stateData, jobKey, reason)
    // Job Non Retryable
    case Event(JobFailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      handleNonRetryableFailure(stateData, jobKey, reason, returnCode)
    // Job Retryable
    case Event(JobFailedRetryableResponse(jobKey, reason, returnCode), _) =>
      handleRetryableFailure(jobKey, reason, returnCode)
    // Aborted? But we're outside of the AbortingState!?? Could happen if
    // - The job was aborted by something external to Cromwell
    // - The job lasted too long (eg PAPI 6 day timeout)
    // Treat it like any other non-retryable failure:
    case Event(JobAbortedResponse(jobKey), stateData) =>
      val cause = new Exception("The job was aborted from outside Cromwell")
      handleNonRetryableFailure(stateData, jobKey, cause)
    // Sub Workflow - sub workflow failures are always non retryable
    case Event(SubWorkflowFailedResponse(jobKey, descendantJobKeys, reason), stateData) =>
      handleNonRetryableFailure(stateData, jobKey, reason, None, descendantJobKeys)
    // Expression evaluation failure
    case Event(ExpressionEvaluationFailedResponse(jobKey, reason), stateData) =>
      handleNonRetryableFailure(stateData, jobKey, reason)

    // Abort command
    case Event(EngineLifecycleActorAbortCommand, data) =>
      handleAbortCommand(data)

    // Misc.
    case Event(RequestValueStore, data) =>
      sender() ! data.valueStore
      stay()
    case Event(EngineStatsActor.JobCountQuery, _) =>
      context.children foreach { _ forward EngineStatsActor.JobCountQuery }
      stay()
    case Event(Terminated(actor), data) =>
      onFailure(actor, new Exception("Actor stopped unexpectedly"))
      stay() using data.removeJobKeyActor(actor)
  }

  onTransition {
    case fromState -> toState if toState.terminal =>
      workflowLogger.debug(s"$tag transitioning from $fromState to $toState. Stopping self.")
      context.stop(self)
    case fromState -> toState =>
      workflowLogger.debug(s"$tag transitioning from $fromState to $toState.")
  }

  /**
    * This is called when a child of this actor throws or terminates. According to Akka, it is safe to access internal data here as long
    * as the supervisor strategy is declared inside the actor (which it is via the StopAndLogSupervisor trait)
    * see https://doc.akka.io/docs/akka/2.5/scala/fault-tolerance.html#creating-a-supervisor-strategy
    */
  override protected def onFailure(actorRef: ActorRef, throwable: => Throwable) = {
    // Both of these Should Never Happen (tm), assuming the state data is set correctly on EJEA creation.
    // If they do, it's a big programmer error and the workflow execution fails.
    val jobKey = stateData.jobKeyActorMappings.getOrElse(actorRef, throw new RuntimeException("Programmer Error: A job or sub workflow actor has terminated but was not assigned a jobKey"))
    val jobStatus = stateData.executionStore.jobStatus(jobKey).getOrElse(throw new RuntimeException(s"Programmer Error: An actor representing ${jobKey.tag} which this workflow is not running has sent up a terminated message."))

    if (!jobStatus.isTerminalOrRetryable) {
      val terminationException = new RuntimeException(s"Unexpected failure or termination of the actor monitoring ${jobKey.tag}", throwable)
      self ! JobFailedNonRetryableResponse(jobKey, terminationException, None)
    }
  }

  private def handleAbortCommand(data: WorkflowExecutionActorData) = {
    workflowLogger.info(s"Aborting workflow")
    // Send the abort to all children
    context.children foreach { _ ! EngineLifecycleActorAbortCommand }
    // As well as all backend singleton actors
    params.backendSingletonCollection.backendSingletonActors.values.flatten.foreach { _ ! BackendSingletonActorAbortWorkflow(workflowIdForLogging) }

    // Only seal the execution store if we're not restarting, otherwise we could miss some jobs that have been started before the
    // restart but are not started yet at this point
    val newData = if (restarting) data else data.sealExecutionStore

    goto(WorkflowExecutionAbortingState) using newData
  }

  private def handleWorkflowSuccessful(data: WorkflowExecutionActorData) = {
    import WorkflowExecutionActor.EnhancedWorkflowOutputs
    import cats.instances.list._
    import cats.syntax.traverse._
    import cromwell.util.JsonFormatting.WomValueJsonFormatter._
    import spray.json._

    def handleSuccessfulWorkflowOutputs(outputs: Map[GraphOutputNode, WomValue]) = {
      val fullyQualifiedOutputs = outputs map {
        case (outputNode, value) => outputNode.identifier.fullyQualifiedName.value -> value
      }
      // Publish fully qualified workflow outputs to log and metadata
      workflowLogger.info(
        s"""Workflow ${workflowDescriptor.callable.name} complete. Final Outputs:
           |${fullyQualifiedOutputs.stripLarge.toJson.prettyPrint}""".stripMargin
      )
      pushWorkflowOutputMetadata(fullyQualifiedOutputs)

      val localOutputs = CallOutputs(outputs map {
        case (outputNode, value) => outputNode.graphOutputPort -> value
      })

      context.parent ! WorkflowExecutionSucceededResponse(data.jobExecutionMap, localOutputs)
      goto(WorkflowExecutionSuccessfulState) using data
    }

    def handleWorkflowOutputsFailure(errors: NonEmptyList[String]) = {
      val exception = new MessageAggregation {
        override def exceptionContext: String = "Workflow output evaluation failed"
        override def errorMessages: Traversable[String] = errors.toList
      }
      context.parent ! WorkflowExecutionFailedResponse(data.jobExecutionMap, exception)
      goto(WorkflowExecutionFailedState)
    }

    // Output ports for the workflow
    val workflowOutputNodes: Set[GraphOutputNode] = workflowDescriptor.callable.graph.outputNodes

    val workflowOutputValuesValidation = workflowOutputNodes
      // Try to find a value for each port in the value store
      .map(outputNode =>
      outputNode -> data.valueStore.get(outputNode.graphOutputPort, None)
    )
      .toList.traverse[IOChecked, (GraphOutputNode, WomValue)]({
      case (name, Some(value)) => value.initialize(data.expressionLanguageFunctions).map(name -> _)
      case (name, None) =>
        s"Cannot find an output value for ${name.identifier.fullyQualifiedName.value}".invalidIOChecked
    })
      // Convert the list of tuples to a Map
      .map(_.toMap)

    workflowOutputValuesValidation
      .map(handleSuccessfulWorkflowOutputs)
      .valueOr(handleWorkflowOutputsFailure)
      .unsafeRunSync()
  }

  def handleWorkflowStalled(data: WorkflowExecutionActorData): State = {
    val notStarted = data.executionStore.unstarted.mkString(System.lineSeparator, System.lineSeparator, "")
    context.parent ! WorkflowExecutionFailedResponse(
      data.jobExecutionMap,
      new Exception(s"Workflow is making no progress but has the following unstarted job keys: $notStarted")
    )
    goto(WorkflowExecutionFailedState)
  }

  private def handleNonRetryableFailure(stateData: WorkflowExecutionActorData,
                                        failedJobKey: JobKey,
                                        reason: Throwable,
                                        returnCode: Option[Int] = None,
                                        jobExecutionMap: JobExecutionMap = Map.empty) = {
    pushFailedCallMetadata(failedJobKey, returnCode, reason, retryableFailure = false)

    val dataWithFailure = stateData.executionFailure(failedJobKey, reason, jobExecutionMap)
    /*
     * If new calls are allowed don't seal the execution store as we want to go as far as possible.
     *
     * Also, if this workflow is restarting, we need to leave all keys in the store to cover for some edge cases.
     * For example, consider the following chain of events:
     *
     * Given 3 calls: A, B and B1, where A and B are independent, and B1 depends on B
     *
     * 1) Start A, Start B
     * 2) B is Successful
     * 3) Start B1
     * 4) Cromwell stops
     * 5) In the meantime, the job for A fails (Cromwell is still down)
     * 5) Cromwell comes back up
     * 6) Recover A, Recover B
     * 7) We get a response saying that A has failed and we end up in this method to handle the failure.
     *
     * At this point if we seal the execution store (remove all NotStarted keys), we'll never get a chance to try to reconnect
     * to B1 which as far as we know might still be running. By leaving the keys in the store, the following can happen:
     *
     * 8) We get a response saying that B was already successful
     * 9) We try to reconnect to B1 (and only reconnect, not recover.
     * We know to do that because we're in Failing state and restarting - see processRunnableJob method)
     * 10) We reconnect to B1 and wait for a terminal status
     * 11) If B1 had a dependency B2, then we'd try to reconnect to it as well, until we find a job that was never started
     * and for which reconnection will fail. When that happens we'll fail the job (see failing state).
     *
     * This guarantees that we'll try to reconnect to all potentially running jobs on restart.
    */
    val newData = if (workflowDescriptor.failureMode.allowNewCallsAfterFailure || restarting) {
      dataWithFailure
    } else {
      dataWithFailure.sealExecutionStore
    }

    // If we're in aborting state, stay there. Even if jobs fail when aborting, the final workflow status should be Aborted.
    val nextState = stateName match {
      case WorkflowExecutionAbortingState => WorkflowExecutionAbortingState
      case _ => WorkflowExecutionFailingState
    }

    goto(nextState) using newData
  }

  private def handleCallAborted(data: WorkflowExecutionActorData, jobKey: JobKey, jobExecutionMap: JobExecutionMap = Map.empty) = {
    pushAbortedCallMetadata(jobKey)
    workflowLogger.info(s"$tag aborted: ${jobKey.tag}")
    val newStateData = data
      .mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Aborted)))
      .addExecutions(jobExecutionMap)

    stay() using newStateData
  }

  private def pushBackendStatusUnknown(jobKey: BackendJobDescriptorKey): Unit = {
    val unknownBackendStatus = MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.BackendStatus), MetadataValue("Unknown"))
    serviceRegistryActor ! PutMetadataAction(unknownBackendStatus)
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey, reason: Throwable, returnCode: Option[Int]) = {
    pushFailedCallMetadata(jobKey, returnCode, reason, retryableFailure = true)

    val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
    workflowLogger.info(s"Retrying job execution for ${newJobKey.tag}")

    // Update current key to RetryableFailure status and add new key with attempt incremented and NotStarted status
    val executionDiff = WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.RetryableFailure, newJobKey -> ExecutionStatus.NotStarted))

    stay() using stateData.mergeExecutionDiff(executionDiff)
  }

  private def handleCallSuccessful(jobKey: JobKey, outputs: CallOutputs, returnCode: Option[Int], data: WorkflowExecutionActorData, jobExecutionMap: JobExecutionMap) = {
    pushSuccessfulCallMetadata(jobKey, returnCode, outputs)
    stay() using data.callExecutionSuccess(jobKey, outputs).addExecutions(jobExecutionMap)
  }

  private def handleDeclarationEvaluationSuccessful(key: ExpressionKey, values: Map[OutputPort, WomValue], data: WorkflowExecutionActorData) = {
    stay() using data.expressionEvaluationSuccess(key, values)
  }

  /**
    * Attempt to start all runnable jobs and return updated state data. This will create a new copy
    * of the state data.
    */
  private def startRunnableNodes(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {
    import keys._

    def updateExecutionStore(diffs: List[WorkflowExecutionDiff], updatedData: WorkflowExecutionActorData): WorkflowExecutionActorData = {
      val notStartedBackendJobs = diffs.flatMap(d => d.executionStoreChanges.collect{
        case (key: BackendJobDescriptorKey, status: ExecutionStatus.NotStarted.type) => (key, status)
      }.keys)
      val notStartedBackendJobsCt = notStartedBackendJobs.size

      // this limits the total max jobs that can be created by a root workflow
      if (data.totalJobsByRootWf.addAndGet(notStartedBackendJobsCt) > TotalMaxJobsPerRootWf) {
        // Since the root workflow tried creating jobs more than the total max jobs allowed per root workflow
        // we fail all the BackendJobDescriptorKey which are in 'Not Started' state, and update the execution
        // store with the status update of remaining keys
        val updatedDiffs = diffs.map(d => d.copy(executionStoreChanges = d.executionStoreChanges -- notStartedBackendJobs))

        notStartedBackendJobs.foreach(key => {
          val errorMsg = s"Job $key failed to be created! Error: Root workflow tried creating ${data.totalJobsByRootWf.get} jobs, which is more than $TotalMaxJobsPerRootWf, the max cumulative jobs allowed per root workflow"
          workflowLogger.error(errorMsg)
          self ! JobFailedNonRetryableResponse(key, new Exception(errorMsg), None)
        })
        updatedData.mergeExecutionDiffs(updatedDiffs)
      }
      else updatedData.mergeExecutionDiffs(diffs)
    }

    val DataStoreUpdate(runnableKeys, statusChanges, updatedData) = data.executionStoreUpdate
    val runnableCalls = runnableKeys.view
      .collect({ case k: BackendJobDescriptorKey => k })
      .groupBy(_.node)
      .map({
        case (node, keys) =>
          val tag = node.fullyQualifiedName
          val shardCount = keys.map(_.index).distinct.size
          if (shardCount == 1) tag
          else s"$tag ($shardCount shards)"
      })
    val mode = if (restarting) "Restarting" else "Starting"
    if (runnableCalls.nonEmpty) workflowLogger.info(s"$mode " + runnableCalls.mkString(", "))

    statusChanges.collect({
      case (jobKey, WaitingForQueueSpace) => pushWaitingForQueueSpaceCallMetadata(jobKey)
    })

    val diffValidation = runnableKeys.traverse[ErrorOr, WorkflowExecutionDiff]({
      case key: BackendJobDescriptorKey => processRunnableJob(key, data)
      case key: SubWorkflowKey => processRunnableSubWorkflow(key, data)
      case key: ConditionalCollectorKey => key.processRunnable(data)
      case key: ConditionalKey => key.processRunnable(data, workflowLogger)
      case key @ ExpressionKey(expr: TaskCallInputExpressionNode, _) => processRunnableTaskCallInputExpression(key, data, expr)
      case key: ExpressionKey => key.processRunnable(data.expressionLanguageFunctions, data.valueStore, self)
      case key: ScatterCollectorKey => key.processRunnable(data)
      case key: ScatteredCallCompletionKey => key.processRunnable(data)
      case key: ScatterKey => key.processRunnable(data, self, MaxScatterWidth)
      case other =>
        workflowLogger.error(s"${other.tag} is not a runnable key")
        WorkflowExecutionDiff.empty.validNel
    })

    // Merge the execution diffs upon success
    diffValidation.map(diffs => updateExecutionStore(diffs, updatedData)).valueOr(errors =>
      throw AggregatedMessageException("Workflow execution failure", errors.toList)
    )
  }

  /*
   * If this ExpressionKey has a TaskCallInputExpressionNode that feeds a task call input, use the backend IoFunctionSet
   * instead of the engine's IoFunctionSet to properly handle `writeFile` or `readFile` invocations.
   */
  private def processRunnableTaskCallInputExpression(key: ExpressionKey,
                                                     data: WorkflowExecutionActorData,
                                                     expressionNode: TaskCallInputExpressionNode): ErrorOr[WorkflowExecutionDiff] = {
    import cats.syntax.either._
    val taskCallNode = expressionNode.taskCallNodeReceivingInput.get(())

    (for {
      backendJobDescriptorKey <- data.executionStore.backendJobDescriptorKeyForNode(taskCallNode) toChecked s"No BackendJobDescriptorKey found for call node $taskCallNode"
      factory <- backendFactoryForTaskCallNode(taskCallNode)
      backendInitializationData = params.initializationData.get(factory.name)
      functions = factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, backendJobDescriptorKey, backendInitializationData, params.ioActor, ioEc)
      diff <- key.processRunnable(functions, data.valueStore, self).toEither
    } yield diff).toValidated
  }

  private def backendFactoryForTaskCallNode(taskCallNode: CommandCallNode): Checked[BackendLifecycleActorFactory] = {
    for {
      name <- workflowDescriptor
        .backendAssignments.get(taskCallNode).toChecked(s"Cannot find an assigned backend for call ${taskCallNode.fullyQualifiedName}")
      factory <- backendFactories.get(name).toChecked(s"Cannot find a backend factory for backend $name")
    } yield factory
  }

  /*
    * Job and Sub Workflow processing
    *
    * Unlike other job keys, those methods are not embedded in the key class itself because they require creating a child actor.
    * While it would be possible to extract those methods from the WEA as well and provide them with an actor factory, the majority of the objects needed to create
    * the children actors are attributes of this class, so it makes more sense to keep the functions here.
   */
  private def processRunnableJob(key: BackendJobDescriptorKey, data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    import cats.syntax.either._
    import common.validation.Checked._

    // Determines the run command to be sent to the BJEA when starting the job
    def runCommand: Checked[BackendJobExecutionActorCommand] = stateName match {
      // In running state, recover if restarting, execute otherwise
      case WorkflowExecutionInProgressState =>
        if (restarting) RecoverJobCommand.validNelCheck else ExecuteJobCommand.validNelCheck

      // In failing state with allowNewCallsAfterFailure, behavior is the same as in running (except we'll fail the workflow at the end)
      case WorkflowExecutionFailingState if workflowDescriptor.failureMode.allowNewCallsAfterFailure =>
        if (restarting) RecoverJobCommand.validNelCheck else ExecuteJobCommand.validNelCheck

      // In failing state without allowNewCallsAfterFailure, try to reconnect if we're restarting
      case WorkflowExecutionFailingState if restarting =>
        ReconnectJobCommand.validNelCheck

      // In aborting state, try to reconnectToAbort if we're restarting
      case WorkflowExecutionAbortingState if restarting =>
        ReconnectToAbortingJobCommand.validNelCheck

      // Any other case is unexpected and is considered a failure (we should not be trying to start a job in any other state)
      case other =>
        s"Cannot start job ${key.tag} in $other state with restarting = $restarting".invalidNelCheck
    }

    (for {
      factory <- backendFactoryForTaskCallNode(key.call)
      command <- runCommand
    } yield startEJEA(key, factory, command)).toValidated
  }

  private def startEJEA(jobKey: BackendJobDescriptorKey,
                        backendLifecycleActorFactory: BackendLifecycleActorFactory,
                        command: BackendJobExecutionActorCommand): WorkflowExecutionDiff = {
    val ejeaName = s"${workflowDescriptor.id}-EngineJobExecutionActor-${jobKey.tag}"
    val backendName = backendLifecycleActorFactory.name
    val backendSingleton = params.backendSingletonCollection.backendSingletonActors(backendName)
    val ejeaProps = EngineJobExecutionActor.props(
      self,
      jobKey,
      workflowDescriptor,
      backendLifecycleActorFactory,
      params.initializationData.get(backendName),
      restarting = params.startState.restarted,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = params.ioActor,
      jobStoreActor = params.jobStoreActor,
      callCacheReadActor = params.callCacheReadActor,
      callCacheWriteActor = params.callCacheWriteActor,
      workflowDockerLookupActor = params.workflowDockerLookupActor,
      jobTokenDispenserActor = params.jobTokenDispenserActor,
      backendSingleton,
      workflowDescriptor.callCachingMode,
      command,
      fileHashCacheActor = params.fileHashCacheActor,
      blacklistCache = params.blacklistCache
    )

    val ejeaRef = context.actorOf(ejeaProps, ejeaName)
    context watch ejeaRef
    pushNewCallMetadata(jobKey, Option(backendName), serviceRegistryActor)
    ejeaRef ! EngineJobExecutionActor.Execute

    WorkflowExecutionDiff(
      executionStoreChanges = Map(jobKey -> ExecutionStatus.QueuedInCromwell),
      jobKeyActorMappings = Map(ejeaRef -> jobKey)
    )
  }

  /*
    * Creates another WEA to process the SubWorkflowKey
   */
  private def processRunnableSubWorkflow(key:SubWorkflowKey, data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    val sweaRef = context.actorOf(
      SubWorkflowExecutionActor.props(
        key,
        workflowDescriptor,
        data.expressionLanguageFunctions,
        backendFactories,
        ioActor = params.ioActor,
        serviceRegistryActor = params.serviceRegistryActor,
        jobStoreActor = params.jobStoreActor,
        subWorkflowStoreActor = params.subWorkflowStoreActor,
        callCacheReadActor = params.callCacheReadActor,
        callCacheWriteActor = params.callCacheWriteActor,
        workflowDockerLookupActor = params.workflowDockerLookupActor,
        jobTokenDispenserActor = params.jobTokenDispenserActor,
        params.backendSingletonCollection,
        params.initializationData,
        params.startState,
        params.rootConfig,
        params.totalJobsByRootWf,
        fileHashCacheActor = params.fileHashCacheActor,
        blacklistCache = params.blacklistCache), s"$workflowIdForLogging-SubWorkflowExecutionActor-${key.tag}"
    )

    context watch sweaRef
    pushNewCallMetadata(key, None, serviceRegistryActor)
    sweaRef ! SubWorkflowExecutionActor.Execute

    WorkflowExecutionDiff(
      executionStoreChanges = Map(key -> ExecutionStatus.QueuedInCromwell),
      jobKeyActorMappings = Map(sweaRef -> key)
    ).validNel
  }
}

object WorkflowExecutionActor {

  /**
    * Rate at which ExecutionHeartBeat events are sent to the the WEA
    */
  val ExecutionHeartBeatInterval = 1 second

  /**
    * States
    */
  sealed trait WorkflowExecutionActorState {
    def terminal = false
  }

  sealed trait WorkflowExecutionActorTerminalState extends WorkflowExecutionActorState {
    override val terminal = true
  }

  case object WorkflowExecutionPendingState extends WorkflowExecutionActorState

  case object WorkflowExecutionInProgressState extends WorkflowExecutionActorState

  case object WorkflowExecutionAbortingState extends WorkflowExecutionActorState

  case object WorkflowExecutionFailingState extends WorkflowExecutionActorState

  case object WorkflowExecutionSuccessfulState extends WorkflowExecutionActorTerminalState

  case object WorkflowExecutionFailedState extends WorkflowExecutionActorTerminalState

  case object WorkflowExecutionAbortedState extends WorkflowExecutionActorTerminalState

  /**
    * Commands
    */
  sealed trait WorkflowExecutionActorCommand

  case object ExecuteWorkflowCommand extends WorkflowExecutionActorCommand

  case object RequestValueStore extends WorkflowExecutionActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowExecutionActorResponse {
    def jobExecutionMap: JobExecutionMap
  }

  case class WorkflowExecutionSucceededResponse(jobExecutionMap: JobExecutionMap, outputs: CallOutputs)
    extends WorkflowExecutionActorResponse {
    override def toString = "WorkflowExecutionSucceededResponse"
  }

  case class WorkflowExecutionAbortedResponse(jobExecutionMap: JobExecutionMap)
    extends WorkflowExecutionActorResponse with EngineLifecycleActorAbortedResponse {
    override def toString = "WorkflowExecutionAbortedResponse"
  }

  final case class WorkflowExecutionFailedResponse(jobExecutionMap: JobExecutionMap, reason: Throwable) extends WorkflowExecutionActorResponse {
    override def toString = "WorkflowExecutionFailedResponse"
  }

  /**
    * Internal control flow messages
    */
  private case class JobInitializationFailed(jobKey: JobKey, throwable: Throwable)

  /**
    * This message triggers the WEA to check for runnable nodes and start them.
    * This operation happens when and only when this message is received.
    * See https://doc.akka.io/docs/akka/2.5/scala/actors.html#timers-scheduled-messages
    */
  private case object ExecutionHeartBeat

  /**
    * Key for the akka timer handling scheduling of the ExecutionHeartBeat
    * See https://doc.akka.io/docs/akka/2.5/scala/actors.html#timers-scheduled-messages
    */
  private case object ExecutionHeartBeatKey

  case class SubWorkflowSucceededResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, outputs: CallOutputs)

  case class SubWorkflowFailedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, reason: Throwable)

  case class SubWorkflowAbortedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap)

  case class WorkflowExecutionActorParams(
                                           workflowDescriptor: EngineWorkflowDescriptor,
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
                                           blacklistCache: Option[BlacklistCache]
                                         )

  def props(workflowDescriptor: EngineWorkflowDescriptor,
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
            blacklistCache: Option[BlacklistCache]): Props = {
    Props(
      WorkflowExecutionActor(
        WorkflowExecutionActorParams(
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
          startState,
          rootConfig,
          totalJobsByRootWf,
          fileHashCacheActor = fileHashCacheActor,
          blacklistCache = blacklistCache
        )
      )
    ).withDispatcher(EngineDispatcher)
  }

  implicit class EnhancedWorkflowOutputs(val outputs: Map[LocallyQualifiedName, WomValue]) extends AnyVal {
    def maxStringLength = 1000

    def stripLarge = outputs map { case (k, v) =>
      val wdlString = v.toWomString

      if (wdlString.length > maxStringLength) (k, WomString(StringUtils.abbreviate(wdlString, maxStringLength)))
      else (k, v)
    }
  }
}
