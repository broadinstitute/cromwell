package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Scope => _, _}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.{AllBackendInitializationData, BackendJobDescriptorKey, JobExecutionMap}
import cromwell.core.Dispatcher._
import cromwell.core.ExecutionIndex._
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.RunnableScopes
import cromwell.engine.workflow.lifecycle.execution.ValueStore.OutputKey
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.keys._
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}
import cromwell.engine.{ContinueWhilePossible, EngineWorkflowDescriptor}
import cromwell.util.StopAndLogSupervisor
import cromwell.webservice.EngineStatsActor
import lenthall.exception.{MessageAggregation, ThrowableAggregation}
import lenthall.util.TryUtil
import lenthall.validation.ErrorOr.ErrorOr
import org.apache.commons.lang3.StringUtils
import wdl._
import wdl.values.WdlArray.WdlArrayLike
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wdl.values.{WdlBoolean, WdlOptionalValue, WdlString, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class WorkflowExecutionActor(workflowDescriptor: EngineWorkflowDescriptor,
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
                                  restarting: Boolean)
  extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] with WorkflowLogging with CallMetadataHelper with StopAndLogSupervisor {

  implicit val ec = context.dispatcher

  override val workflowIdForLogging = workflowDescriptor.id
  override val workflowIdForCallMetadata = workflowDescriptor.id

  private val tag = s"WorkflowExecutionActor [UUID(${workflowDescriptor.id.shortString})]"

  private var checkRunnableCancellable: Option[Cancellable] = None

  private val backendFactories = TryUtil.sequenceMap(workflowDescriptor.backendAssignments.values.toSet[String] map { backendName =>
    backendName -> CromwellBackends.backendLifecycleFactoryActorByName(backendName)
  } toMap) recover {
    case e => throw new RuntimeException("Could not instantiate backend factories", e)
  } get

  startWith(
    WorkflowExecutionPendingState,
    WorkflowExecutionActorData(
      workflowDescriptor,
      executionStore = ExecutionStore(workflowDescriptor.backendDescriptor.workflow),
      backendJobExecutionActors = Map.empty,
      engineCallExecutionActors = Map.empty,
      subWorkflowExecutionActors = Map.empty,
      downstreamExecutionMap = Map.empty,
      valueStore = ValueStore.initialize(workflowDescriptor.knownValues)
    )
  )

  when(WorkflowExecutionPendingState) {
    case Event(ExecuteWorkflowCommand, _) =>
      // TODO WOM: Remove this when conditional and sub workflows are supported. It prevents the workflow from hanging
      if(workflowDescriptor.namespace.innerGraph.nodes.collectFirst({
        case  _: WorkflowCallNode => true
      }).nonEmpty) {
        context.parent ! WorkflowExecutionFailedResponse(Map.empty, new Exception("Scatter and Conditional not supported yet"))
        goto(WorkflowExecutionFailedState)
      }
      scheduleStartRunnableCalls()
      goto(WorkflowExecutionInProgressState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(RequestValueStore, data) =>
      sender() ! data.valueStore
      stay()
    case Event(CheckRunnable, data) => handleCheckRunnable(data)

    case Event(JobStarting(jobKey), stateData) =>
      pushStartingCallMetadata(jobKey)
      stay() using stateData
        .mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Starting)))
    case Event(JobRunning(key, inputs, callExecutionActor), stateData) =>
      pushRunningCallMetadata(key, inputs)
      stay() using stateData
        .addCallExecutionActor(key, callExecutionActor)
        .mergeExecutionDiff(WorkflowExecutionDiff(Map(key -> ExecutionStatus.Running)))

    //Success
    // Job
    case Event(r: JobSucceededResponse, stateData) =>
      pushSuccessfulCallMetadata(r.jobKey, r.returnCode, r.jobOutputs)
      handleCallSuccessful(r.jobKey, r.jobOutputs, stateData, Map.empty)
    // Sub Workflow
    case Event(SubWorkflowSucceededResponse(jobKey, descendantJobKeys, callOutputs), stateData) =>
      pushSuccessfulCallMetadata(jobKey, None, callOutputs)
      handleCallSuccessful(jobKey, callOutputs, stateData, descendantJobKeys)
    // Expression
    case Event(ExpressionEvaluationSucceededResponse(jobKey, callOutputs), stateData) =>
      handleDeclarationEvaluationSuccessful(jobKey, callOutputs, stateData)
    // Conditional
    case Event(BypassedCallResults(callOutputs), stateData) =>
      handleCallBypassed(callOutputs, stateData)
    case Event(BypassedDeclaration(declKey), stateData) =>
      handleDeclarationEvaluationSuccessful(declKey, WdlOptionalValue.none(declKey.womType), stateData)

    // Failure
    // Initialization
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      pushFailedCallMetadata(jobKey, None, reason, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, reason, Map.empty)
    // Job Non Retryable
    case Event(JobFailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedCallMetadata(jobKey, returnCode, reason, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, reason, Map.empty)
    // Aborted? But we're outside of the AbortingState!?? Could happen if
    // - The job was aborted by something external to Cromwell
    // - The job lasted too long (eg JES 6 day timeout)
    // - We've reconnected to an aborting job (some sort of shutdown race condition?)
    // Treat it like any other non-retryable failure:
    case Event(AbortedResponse(jobKey), stateData) =>
      val cause = new Exception("The job was aborted from outside Cromwell")
      pushFailedCallMetadata(jobKey, None, cause, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, cause, Map.empty)
    // Job Retryable
    case Event(JobFailedRetryableResponse(jobKey, reason, returnCode), _) =>
      pushFailedCallMetadata(jobKey, None, reason, retryableFailure = true)
      handleRetryableFailure(jobKey, reason, returnCode)
    // Sub Workflow - sub workflow failures are always non retryable
    case Event(SubWorkflowFailedResponse(jobKey, descendantJobKeys, reason), stateData) =>
      pushFailedCallMetadata(jobKey, None, reason, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, reason, descendantJobKeys)
    case Event(ExpressionEvaluationFailedResponse(jobKey, reason), stateData) =>
      handleDeclarationEvaluationFailure(jobKey, reason, stateData)
  }

  when(WorkflowExecutionAbortingState) {
    case Event(AbortedResponse(jobKey), stateData) =>
      pushAbortedCallMetadata(jobKey)
      handleCallAborted(stateData, jobKey, Map.empty)
    case Event(SubWorkflowAbortedResponse(jobKey, executedKeys), stateData) =>
      pushAbortedCallMetadata(jobKey)
      handleCallAborted(stateData, jobKey, executedKeys)
    case Event(SubWorkflowSucceededResponse(subKey, executedKeys, _), stateData) =>
      pushAbortedCallMetadata(subKey)
      handleCallAborted(stateData, subKey, executedKeys)
    case Event(r: JobSucceededResponse, stateData) =>
      pushAbortedCallMetadata(r.jobKey)
      handleCallAborted(stateData, r.jobKey, Map.empty)
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

  private def scheduleStartRunnableCalls() = {
    checkRunnableCancellable = Option(context.system.scheduler.scheduleOnce(SweepInterval, self, CheckRunnable))
  }

  override def postStop(): Unit = {
    checkRunnableCancellable foreach { _.cancel() }
    super.postStop()
  }

  private def handleTerminated(actorRef: ActorRef) = {
    // Both of these Should Never Happen (tm), assuming the state data is set correctly on EJEA creation.
    // If they do, it's a big programmer error and the workflow execution fails.
    val jobKey = stateData.engineCallExecutionActors.getOrElse(actorRef, throw new RuntimeException("Programmer Error: An EJEA has terminated but was not assigned a jobKey"))
    val jobStatus = stateData.executionStore.jobStatus(jobKey).getOrElse(throw new RuntimeException("Programmer Error: An EJEA representing a jobKey which this workflow is not running has sent up a terminated message."))

    if (!jobStatus.isTerminalOrRetryable) {
      val terminationException = getFailureCause(actorRef) match {
        case Some(e) => new RuntimeException("Unexpected failure (or early exit) in EJEA.", e)
        case None => new RuntimeException(s"Unexpected failure (or early exit) in EJEA $actorRef (root cause not captured).")
      }
      self ! JobFailedNonRetryableResponse(jobKey, terminationException, None)
    }

    stay
  }

  whenUnhandled {
    case Event(CheckRunnable, _) => stay()
    case Event(Terminated(actorRef), stateData) => handleTerminated(actorRef) using stateData.removeEngineJobExecutionActor(actorRef)
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      if (stateData.hasRunningActors) {
        log.info(s"$tag: Abort received. " +
          s"Aborting ${stateData.backendJobExecutionActors.size} Job Execution Actors" +
          s" and ${stateData.subWorkflowExecutionActors.size} Sub Workflow Execution Actors"
        )
        stateData.backendJobExecutionActors.values foreach { _ ! AbortJobCommand }
        stateData.subWorkflowExecutionActors.values foreach { _ ! EngineLifecycleActorAbortCommand }
        goto(WorkflowExecutionAbortingState)
      } else {
        goto(WorkflowExecutionAbortedState)
      }
    case Event(EngineStatsActor.JobCountQuery, data) =>
      sender ! EngineStatsActor.JobCount(data.backendJobExecutionActors.size)
      data.subWorkflowExecutionActors.values foreach { _ forward EngineStatsActor.JobCountQuery }
      stay()
    case unhandledMessage =>
      workflowLogger.warn(s"$tag received an unhandled message: ${unhandledMessage.event} in state: $stateName")
      stay()
  }

  onTransition {
    case fromState -> toState if toState.terminal =>
      workflowLogger.debug(s"$tag transitioning from $fromState to $toState. Stopping self.")
      context.stop(self)
    case fromState -> toState =>
      workflowLogger.debug(s"$tag transitioning from $fromState to $toState.")
  }

  onTransition {
    case _ -> WorkflowExecutionAbortedState =>
      context.parent ! WorkflowExecutionAbortedResponse(nextStateData.jobExecutionMap)
  }

  private def handleNonRetryableFailure(stateData: WorkflowExecutionActorData, failedJobKey: JobKey, reason: Throwable, jobExecutionMap: JobExecutionMap) = {
    val newData = stateData
      .removeCallExecutionActor(failedJobKey)
      .addExecutions(jobExecutionMap)

    handleExecutionFailure(failedJobKey, newData, reason, jobExecutionMap)
  }

  private def handleDeclarationEvaluationFailure(declarationKey: ExpressionKey, reason: Throwable, stateData: WorkflowExecutionActorData) = {
    handleExecutionFailure(declarationKey, stateData, reason, Map.empty)
  }

  private def handleExecutionFailure(failedJobKey: JobKey, data: WorkflowExecutionActorData, reason: Throwable, jobExecutionMap: JobExecutionMap) = {
    val newData = data.executionFailed(failedJobKey)

    if (workflowDescriptor.failureMode == ContinueWhilePossible) {
      newData.workflowCompletionStatus match {
        case Some(completionStatus) if completionStatus == Failed =>
          context.parent ! WorkflowExecutionFailedResponse(newData.jobExecutionMap, reason)
          goto(WorkflowExecutionFailedState) using newData
        case _ =>
          stay() using newData
      }
    } else {
      context.parent ! WorkflowExecutionFailedResponse(newData.jobExecutionMap, reason)
      goto(WorkflowExecutionFailedState) using newData
    }
  }

  private def handleWorkflowSuccessful(data: WorkflowExecutionActorData) = {
    import WorkflowExecutionActor.EnhancedWorkflowOutputs
    import cats.instances.list._
    import cats.syntax.traverse._
    import cromwell.util.JsonFormatting.WdlValueJsonFormatter._
    import spray.json._

    def handleSuccessfulWorkflowOutputs(outputs: Map[WomIdentifier, WdlValue]) = {
      val fullyQualifiedOutputs = outputs map {
        case (identifier, value) => identifier.fullyQualifiedName.value -> value
      }
      // Publish fully qualified workflow outputs to log and metadata
      workflowLogger.info(
        s"""Workflow ${workflowDescriptor.workflow.name} complete. Final Outputs:
           |${fullyQualifiedOutputs.stripLarge.toJson.prettyPrint}""".stripMargin
      )
      pushWorkflowOutputMetadata(fullyQualifiedOutputs)

      // Use local names so they can be used in outer workfows if this is a sub workflow
      val localOutputs = outputs map {
        case (identifier, value) => identifier.localName.value -> JobOutput(value)
      }

      context.parent ! WorkflowExecutionSucceededResponse(data.jobExecutionMap, localOutputs)
      goto(WorkflowExecutionSuccessfulState) using data
    }

    workflowDescriptor.namespace.innerGraph.outputNodes
      .flatMap(_.outputPorts)
      .map(op => op.identifier -> data.valueStore.get(op, None)).toList
      .traverse[ErrorOr, (WomIdentifier, WdlValue)]({
      case (name, Some(value)) => (name -> value).validNel
      case (name, None) => s"Cannot find an output value for ${name.fullyQualifiedName.value}".invalidNel
    }).map(validOutputs => handleSuccessfulWorkflowOutputs(validOutputs.toMap))
      .valueOr { errors =>
        val exception = new MessageAggregation {
          override def exceptionContext: String = "Workflow output evaluation failed"
          override def errorMessages: Traversable[String] = errors.toList
        }
        context.parent ! WorkflowExecutionFailedResponse(data.jobExecutionMap, exception)
        goto(WorkflowExecutionFailedState)
      }
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey, reason: Throwable, returnCode: Option[Int]) = {
    val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
    workflowLogger.info(s"Retrying job execution for ${newJobKey.tag}")
    /*  Currently, we update the status of the old key to RetryableFailure, and add a new entry (with the #attempts incremented by 1)
      * to the execution store with status as NotStarted. This allows startRunnableCalls to re-execute this job */
    val executionDiff = WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.RetryableFailure, newJobKey -> ExecutionStatus.NotStarted))
    val newData = stateData.mergeExecutionDiff(executionDiff).removeCallExecutionActor(jobKey)
    stay() using newData
  }

  private def handleCallSuccessful(jobKey: JobKey, outputs: CallOutputs, data: WorkflowExecutionActorData, jobExecutionMap: JobExecutionMap) = {
    stay() using data.callExecutionSuccess(jobKey, outputs).addExecutions(jobExecutionMap)
  }

  private def handleDeclarationEvaluationSuccessful(key: ExpressionKey, value: WdlValue, data: WorkflowExecutionActorData) = {
    stay() using data.expressionEvaluationSuccess(key, value)
  }

  private def handleCallBypassed(callOutputs: Map[CallKey, CallOutputs], data: WorkflowExecutionActorData) = {
    def foldFunc(d: WorkflowExecutionActorData, output: (CallKey, CallOutputs)) = d.callExecutionSuccess(output._1, output._2)

    val updatedData = callOutputs.foldLeft(data)(foldFunc)
    stay() using updatedData
  }

  private def handleCheckRunnable(data: WorkflowExecutionActorData) = {
    data.workflowCompletionStatus match {
      case Some(ExecutionStatus.Done) =>
        handleWorkflowSuccessful(data)
      case Some(_) =>
        context.parent ! WorkflowExecutionFailedResponse(data.jobExecutionMap, new Exception("One or more jobs failed in ContinueWhilePossible mode"))
        goto(WorkflowExecutionFailedState) using data
      case _ =>
        scheduleStartRunnableCalls()
        if (data.hasNewRunnables) stay() using startRunnableScopes(data) else stay()
    }
  }

  private def handleCallAborted(data: WorkflowExecutionActorData, jobKey: JobKey, jobExecutionMap: JobExecutionMap) = {
    workflowLogger.info(s"$tag job aborted: ${jobKey.tag}")
    val newStateData = data
      .mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Aborted)))
      .removeCallExecutionActor(jobKey)
      .addExecutions(jobExecutionMap)
    if (!newStateData.hasRunningActors) {
      workflowLogger.info(s"$tag all jobs aborted")
      goto(WorkflowExecutionAbortedState)
    } else {
      stay() using newStateData
    }
  }

  /**
    * Attempt to start all runnable jobs and return updated state data.  This will create a new copy
    * of the state data.
    */
  private def startRunnableScopes(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {
    val RunnableScopes(runnableScopes, truncated) = data.executionStore.runnableScopes
    val runnableCalls = runnableScopes.view collect { case k if k.node.isInstanceOf[CallNode] => k } sortBy { k =>
      (k.node.fullyQualifiedName, k.index.getOrElse(-1)) } map { _.tag }

    if (runnableCalls.nonEmpty) workflowLogger.info("Starting calls: " + runnableCalls.mkString(", "))

    // Each process returns a Try[WorkflowExecutionDiff], which, upon success, contains potential changes to be made to the execution store.
    val diffs = runnableScopes map { scope =>
      scope -> Try(scope match {
        case k: CallKey if isInBypassedScope(k, data) => processBypassedScope(k, data)
        case k: ExpressionKey if isInBypassedScope(k, data) => processBypassedScope(k, data)
        case k: BackendJobDescriptorKey => processRunnableJob(k, data)
        case k: ScatterKey => processRunnableScatter(k, data, isInBypassedScope(k, data))
        case k: ConditionalKey => processRunnableConditional(k, data)
        case k: ScatterCollectorKey => processRunnableScatterCollector(k, data, isInBypassedScope(k, data))
        case k: ConditionalCollectorKey => processRunnableConditionalCollector(k, data, isInBypassedScope(k, data))
        case k: SubWorkflowKey => processRunnableSubWorkflow(k, data)
        case k: ExpressionKey => processRunnableExpression(k, data)
        case k => Failure(new UnsupportedOperationException(s"Unknown entry in execution store: ${k.tag}"))
      }).flatten
    } map {
      case (_, Success(value)) => Success(value)
      case (scope, Failure(throwable)) =>
        self ! JobInitializationFailed(scope, throwable)
        Failure(throwable)
    } collect {
      /*
      NOTE: This is filtering out all errors and only returning the successes, but only after the map above sent a
      message that something is wrong.

      We used to throw an aggregation exception of all the collected errors, but _nothing_ in cromwell is actually
      expecting that. Thus the workflows were being left in a Running state, jobs were left dispatched, etc.

      Meanwhile this actor and its children died or restarted, and one couldn't even attempt to abort the other jobs.

      Now, in the previous map, we send a message to ourselves about _every_ failure. But this method does not attempt
      to further process the errors. The failures enqueue in the actor mailbox, and are handled by this actor's receive.

      At the moment, there is an issue in how this actor handles failure messages. That issue is tracked in:
      https://github.com/broadinstitute/cromwell/issues/2029

      Separately, we may also want to institute better supervision of actors, in general. But just throwing an exception
      here doesn't actually force the correct handling.

      See also:
      https://github.com/broadinstitute/cromwell/issues/1414
      https://github.com/broadinstitute/cromwell/issues/1874
      */
      case Success(value) => value
    }

    // Update the metadata for the jobs we just sent to EJEAs (they'll start off queued up waiting for tokens):
    pushQueuedCallMetadata(diffs)
    val newData = data.mergeExecutionDiffs(diffs)
    if (truncated || diffs.exists(_.containsNewEntry)) newData else newData.resetCheckRunnable
  }

  private def isInBypassedScope(jobKey: JobKey, data: WorkflowExecutionActorData) = {
    data.executionStore.isInBypassedConditional(jobKey)
  }

  private def processBypassedScope(jobKey: JobKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    Success(
      WorkflowExecutionDiff(
        executionStoreChanges = Map(jobKey -> ExecutionStatus.Bypassed),
        valueStoreAdditions = bypassedScopeResults(jobKey)
      )
    )
  }

  private def bypassedScopeResults(jobKey: JobKey): Map[OutputKey, WdlOptionalValue] = {
    jobKey.node.outputPorts.map({ outputPort =>
      OutputKey(outputPort, jobKey.index) ->  WdlOptionalValue.none(outputPort.womType)
    }).toMap
  }

  private def processRunnableExpression(expression: ExpressionKey, data: WorkflowExecutionActorData) = {
    import lenthall.validation.ErrorOr._

    expression.upstreamPorts.traverseValues(resolve(expression, data)) map { lookup =>
      // Send a message to self in case we decide to change evaluate to return asynchronously, if we don't we could
      // directly add the value to the value store in the execution diff
      expression.node.evaluateAndCoerce(lookup, data.expressionLanguageFunctions) match {
        case Right(result) => self ! ExpressionEvaluationSucceededResponse(expression, result)
        case Left(f) => self ! ExpressionEvaluationFailedResponse(expression, new RuntimeException(f.toList.mkString(", ")))
      }
    } valueOr { f =>
      self ! ExpressionEvaluationFailedResponse(expression, new RuntimeException(f.toList.mkString(", ")))
    }

    Success(WorkflowExecutionDiff(Map(expression -> ExecutionStatus.Running)))
  }

  /**
    * Attempts to resolve an output port to a known value.
    * Curried for convenience.
    */
  private def resolve(jobKey: JobKey, data: WorkflowExecutionActorData)(outputPort: OutputPort): ErrorOr[WdlValue] = {
    // If the node this port belongs to is a ScatterVariableNode then we want the item at the right index in the array
    def forScatterVariable: ErrorOr[WdlValue] = data.valueStore.get(outputPort, None) match {
      // Try to find the element at "jobIndex" in the array value stored for the outputPort, any other case is a failure
      case Some(wdlValue: WdlArrayLike) =>
        jobKey.index match {
          case Some(jobIndex) =>
            wdlValue.asArray.value.lift(jobIndex)
              .map(_.validNel)
              .getOrElse(s"Shard index $jobIndex exceeds scatter array length: ${wdlValue.asArray.value.size}".invalidNel)
          case None => s"Unsharded execution key ${jobKey.tag} references a scatter variable: ${outputPort.identifier.fullyQualifiedName}".invalidNel
        }
      case Some(other) => s"Value for scatter collection ${outputPort.identifier.fullyQualifiedName} is not an array: ${other.wdlType.toWdlString}".invalidNel
      case None => s"Can't find a value for scatter collection ${outputPort.identifier.fullyQualifiedName}".invalidNel
    }
    
    // Just look it up in the store
    def forNormalNode(index: ExecutionIndex) = data.valueStore.get(outputPort, index) match {
      case Some(value) => value.validNel
      case None => s"Can't find a value for ${outputPort.name}".invalidNel
    }

    outputPort.graphNode match {
      case _: ScatterVariableNode => forScatterVariable
        // OuterGraphInputNodes are not indexed (even if this jobKey is a shard, the node is outside the scatter)
      case _: OuterGraphInputNode => forNormalNode(None)
      case _ => forNormalNode(jobKey.index)
    }
  }

  private def processRunnableJob(jobKey: BackendJobDescriptorKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        val exception = new IllegalStateException(s"$tag $message")
        workflowLogger.error(exception, s"$tag $message")
        throw exception
      case Some(backendName) =>
        backendFactories.get(backendName) match {
          case Some(factory) =>
            val ejeaName = s"${workflowDescriptor.id}-EngineJobExecutionActor-${jobKey.tag}"
            val backendSingleton = backendSingletonCollection.backendSingletonActors(backendName)
            val ejeaProps = EngineJobExecutionActor.props(
              self, jobKey, workflowDescriptor, factory, initializationData.get(backendName), restarting,
              serviceRegistryActor = serviceRegistryActor,
              ioActor = ioActor,
              jobStoreActor = jobStoreActor,
              callCacheReadActor = callCacheReadActor,
              callCacheWriteActor = callCacheWriteActor,
              workflowDockerLookupActor = workflowDockerLookupActor,
              jobTokenDispenserActor = jobTokenDispenserActor,
              backendSingleton, backendName, workflowDescriptor.callCachingMode)
            val ejeaRef = context.actorOf(ejeaProps, ejeaName)
            context watch ejeaRef
            pushNewCallMetadata(jobKey, Option(backendName))
            ejeaRef ! EngineJobExecutionActor.Execute
            Success(WorkflowExecutionDiff(
              executionStoreChanges = Map(jobKey -> ExecutionStatus.QueuedInCromwell),
              engineJobExecutionActorAdditions = Map(ejeaRef -> jobKey)))
          case None =>
            throw WorkflowExecutionException(NonEmptyList.of(new Exception(s"Could not get BackendLifecycleActor for backend $backendName")))
        }
    }
  }

  private def processRunnableSubWorkflow(key: SubWorkflowKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    val sweaRef = context.actorOf(
      SubWorkflowExecutionActor.props(key,
        data.workflowDescriptor,
        data.expressionLanguageFunctions,
        backendFactories,
        ioActor = ioActor,
        serviceRegistryActor = serviceRegistryActor,
        jobStoreActor = jobStoreActor,
        subWorkflowStoreActor = subWorkflowStoreActor,
        callCacheReadActor = callCacheReadActor,
        callCacheWriteActor = callCacheWriteActor,
        workflowDockerLookupActor = workflowDockerLookupActor,
        jobTokenDispenserActor = jobTokenDispenserActor,
        backendSingletonCollection, initializationData, restarting), s"SubWorkflowExecutionActor-${key.tag}"
    )

    context watch sweaRef
    pushNewCallMetadata(key, None)
    sweaRef ! SubWorkflowExecutionActor.Execute

    Success(WorkflowExecutionDiff(executionStoreChanges = Map(key -> ExecutionStatus.QueuedInCromwell),
      engineJobExecutionActorAdditions = Map(sweaRef -> key)))
  }

  private def processRunnableConditional(conditionalKey: ConditionalKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    val conditionOutputPort = conditionalKey.node.conditionExpression.singleExpressionOutputPort

    data.valueStore.get(conditionOutputPort, conditionalKey.index) match {
      case Some(b: WdlBoolean) =>
        val conditionalStatus = if (b.value) ExecutionStatus.Done else ExecutionStatus.Bypassed
        Success(WorkflowExecutionDiff(conditionalKey.populate + (conditionalKey -> conditionalStatus)))
      case Some(v: WdlValue) => Failure(new RuntimeException(
        s"'if' condition must evaluate to a boolean but instead got ${v.wdlType.toWdlString}"))
      case None => Failure(
        new RuntimeException(s"Could not find the boolean value for conditional ${conditionalKey.tag}")
      )
    }
  }

  private def processRunnableScatter(scatterKey: ScatterKey, data: WorkflowExecutionActorData, bypassed: Boolean): Try[WorkflowExecutionDiff] = {
    if (bypassed) {
      Success(WorkflowExecutionDiff(scatterKey.populate(0) + (scatterKey -> ExecutionStatus.Bypassed)))
    } else {
      val collectionOutputPort = scatterKey.node.scatterCollectionExpressionNode.singleExpressionOutputPort

      data.valueStore.get(collectionOutputPort, None) map {
        case WdlArrayLike(arrayLike) =>
          Success(
            WorkflowExecutionDiff(
              // Add the new shards + collectors, and set the scatterVariable and scatterKey to Done
              executionStoreChanges = scatterKey.populate(arrayLike.value.size) ++ Map (
                ScatterVariableInputKey(scatterKey.node.scatterVariableInnerGraphInputNode, arrayLike) -> Done,
                scatterKey -> Done
              ),
              // Add scatter variable arrayLike to the value store
              valueStoreAdditions = Map (
                OutputKey(scatterKey.node.scatterVariableInnerGraphInputNode.singleOutputPort, None) -> arrayLike
              )
            )
          )
        case v: WdlValue =>
          Failure(new RuntimeException(s"Scatter collection must evaluate to an array but instead got ${v.wdlType.toWdlString}"))
      } getOrElse {
        Failure(new RuntimeException(s"Could not find an array value for scatter ${scatterKey.tag}"))
      }
    }
  }

  /**
    * Collects all shards and add them to the value store
    */
  private def processRunnableScatterCollector(collector: ScatterCollectorKey, data: WorkflowExecutionActorData, isInBypassed: Boolean): Try[WorkflowExecutionDiff] = {
    data.valueStore.collectShards(collector) match {
      case Invalid(e) => Failure(new RuntimeException with MessageAggregation {
        override def exceptionContext: String = s"Failed to collect output shards for node ${collector.tag}"
        override def errorMessages: Traversable[String] = e.toList
      })
      case Valid(outputs) =>
        // TODO WOM: fix
        //        val adjustedOutputs: CallOutputs = if (isInBypassed) {
        //          outputs map {
        //            case (outputKey, value) => outputKey -> WdlOptionalValue.none(output._2.wdlValue.wdlType
        //          }
        //        } else outputs
        Success(WorkflowExecutionDiff(
          executionStoreChanges = Map(collector -> ExecutionStatus.Done),
          valueStoreAdditions = outputs
        ))
    }
  }

  private def processRunnableConditionalCollector(collector: ConditionalCollectorKey, data: WorkflowExecutionActorData, isInBypassed: Boolean): Try[WorkflowExecutionDiff] = {
    data.valueStore.collectConditional(collector, isInBypassed) match {
      case Invalid(e) => Failure(new RuntimeException with MessageAggregation {
        override def exceptionContext: String = s"Failed to collect conditional output value for node ${collector.tag}"
        override def errorMessages: Traversable[String] = e.toList
      })
      case Valid(outputs) => Success(WorkflowExecutionDiff(
        executionStoreChanges = Map(collector -> ExecutionStatus.Done),
        valueStoreAdditions = outputs
      ))
    }
  }
}

object WorkflowExecutionActor {

  val SweepInterval = 1 second

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

  private case class ExpressionEvaluationSucceededResponse(declarationKey: ExpressionKey, value: WdlValue)

  private case object CheckRunnable

  private[execution] sealed trait BypassedScopeResults

  private case class BypassedCallResults(callOutputs: Map[CallKey, CallOutputs]) extends BypassedScopeResults
  private case class BypassedDeclaration(declaration: ExpressionKey) extends BypassedScopeResults

  private case class ExpressionEvaluationFailedResponse(declarationKey: ExpressionKey, reason: Throwable)

  case class SubWorkflowSucceededResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, outputs: CallOutputs)

  case class SubWorkflowFailedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, reason: Throwable)

  case class SubWorkflowAbortedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap)

  case class WorkflowExecutionException[T <: Throwable](exceptions: NonEmptyList[T]) extends ThrowableAggregation {
    override val throwables = exceptions.toList
    override val exceptionContext = s"WorkflowExecutionActor"
  }

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
            restarting: Boolean): Props = {
    Props(WorkflowExecutionActor(workflowDescriptor,
      ioActor = ioActor,
      serviceRegistryActor = serviceRegistryActor,
      jobStoreActor = jobStoreActor,
      subWorkflowStoreActor = subWorkflowStoreActor,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      workflowDockerLookupActor = workflowDockerLookupActor,
      jobTokenDispenserActor = jobTokenDispenserActor,
      backendSingletonCollection, initializationData, restarting)).withDispatcher(EngineDispatcher)
  }

  implicit class EnhancedWorkflowOutputs(val outputs: Map[LocallyQualifiedName, WdlValue]) extends AnyVal {
    def maxStringLength = 1000

    def stripLarge = outputs map { case (k, v) =>
      val wdlString = v.toWdlString

      if (wdlString.length > maxStringLength) (k, WdlString(StringUtils.abbreviate(wdlString, maxStringLength)))
      else (k, v)
    }
  }
}
