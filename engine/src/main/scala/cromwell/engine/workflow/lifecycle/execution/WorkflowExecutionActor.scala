package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Scope => _, _}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.{AllBackendInitializationData, BackendJobDescriptorKey, JobExecutionMap}
import cromwell.core.CromwellGraphNode._
import cromwell.core.Dispatcher._
import cromwell.core.ExecutionIndex._
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.lifecycle.execution.ExecutionStore.RunnableScopes
import cromwell.engine.workflow.lifecycle.execution.OutputStore.OutputKey
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor._
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}
import cromwell.engine.{ContinueWhilePossible, EngineWorkflowDescriptor}
import cromwell.util.StopAndLogSupervisor
import cromwell.webservice.EngineStatsActor
import lenthall.exception.ThrowableAggregation
import lenthall.util.TryUtil
import lenthall.validation.ErrorOr.ErrorOr
import org.apache.commons.lang3.StringUtils
import wdl._
import wdl.values.{WdlOptionalValue, WdlString, WdlValue}
import wom.expression.IoFunctionSet
import wom.graph.GraphNodePort.{GraphNodeOutputPort, InputPort, OutputPort}
import wom.graph._

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
      executionStore = ExecutionStore(workflowDescriptor.backendDescriptor.workflow, workflowDescriptor.knownValues),
      backendJobExecutionActors = Map.empty,
      engineCallExecutionActors = Map.empty,
      subWorkflowExecutionActors = Map.empty,
      downstreamExecutionMap = Map.empty,
      outputStore = OutputStore.initialize(workflowDescriptor.knownValues)
    )
  )

  when(WorkflowExecutionPendingState) {
    case Event(ExecuteWorkflowCommand, _) =>
      scheduleStartRunnableCalls()
      goto(WorkflowExecutionInProgressState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(RequestOutputStore, data) =>
      sender() ! data.outputStore
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
    // Scatter
    case Event(ScatterCollectionSucceededResponse(jobKey, callOutputs), stateData) =>
      handleCallSuccessful(jobKey, callOutputs, stateData, Map.empty)
    // Declaration
    case Event(DeclarationEvaluationSucceededResponse(jobKey, callOutputs), stateData) =>
      handleDeclarationEvaluationSuccessful(jobKey, callOutputs, stateData)
    // Conditional
    case Event(BypassedCallResults(callOutputs), stateData) =>
      handleCallBypassed(callOutputs, stateData)
    case Event(BypassedDeclaration(declKey), stateData) =>
      handleDeclarationEvaluationSuccessful(declKey, WdlOptionalValue.none(declKey.node.womType), stateData)

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
    case Event(DeclarationEvaluationFailedResponse(jobKey, reason), stateData) =>
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

  override def postStop() = {
    checkRunnableCancellable foreach { _.cancel() }
    super.postStop()
  }

  def handleTerminated(actorRef: ActorRef) = {
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
    import cromwell.util.JsonFormatting.WdlValueJsonFormatter._
    import spray.json._

    case class ResponseAndFinalState(response: WorkflowExecutionActorResponse,
                                     finalState: WorkflowExecutionActorTerminalState)

    def shouldFilterOutputPort(outputPort: OutputPort): Boolean = {
      outputPort match {
        case gnop: GraphNodeOutputPort =>
          gnop.graphNode match {
            case _: RequiredGraphInputNode => true // Yes RequiredGraphInputNodes have output ports, but we don't want
                                                   // to recycle them back to outputs.
            case _: ExpressionNode => true // Expressions whose output ports feed into other nodes.
            case _ => false // Default to not filtering.
          }
        case _ => false // Default to not filtering.
      }
    }

    // TODO WOM: workflow outputs ? Just dump the outputStore for now...
    // For logging and metadata
    val outputs: Map[String, WdlValue] = data.outputStore.store collect {
      case (OutputKey(outputPort, _), value) if !shouldFilterOutputPort(outputPort) => s"${outputPort.fullyQualifiedName}" -> value
    }

    val workflowScopeOutputs: Map[String, WdlValue] = outputs map {
      case (key, value) => s"${workflowDescriptor.workflow.name}.$key" -> value
    }

    workflowLogger.info(
      s"""Workflow ${workflowDescriptor.workflow.name} complete. Final Outputs:
         |${workflowScopeOutputs.stripLarge.toJson.prettyPrint}""".stripMargin
    )
    pushWorkflowOutputMetadata(workflowScopeOutputs)

    // For cromwell internal storage of outputs
    val unqualifiedWorkflowOutputs = outputs map {
      // JobOutput is poorly named here - a WorkflowOutput type would be better
      case (output, value) => output -> JobOutput(value)
    }
    val responseAndState = ResponseAndFinalState(
      WorkflowExecutionSucceededResponse(data.jobExecutionMap, unqualifiedWorkflowOutputs),
      WorkflowExecutionSuccessfulState)


    context.parent ! responseAndState.response
    goto(responseAndState.finalState) using data
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
    stay() using data.declarationEvaluationSuccess(key, value)
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
        case k: CollectorKey => processRunnableCollector(k, data, isInBypassedScope(k, data))
        case k: SubWorkflowKey => processRunnableSubWorkflow(k, data)
        case k: DynamicDeclarationKey => processRunnableDynamicDeclaration(k, data)
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
    //    val result = jobKey.scope.ancestry.exists {
    //      case i: If => data.executionStore.isBypassedConditional(jobKey, i)
    //      case _ => false
    //    }
    //    result
    false
  }

  def processBypassedScope(jobKey: JobKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    self ! bypassedScopeResults(jobKey)
    Success(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Running)))
  }

  def bypassedScopeResults(jobKey: JobKey): BypassedScopeResults = {
    //    jobKey match {
    //      case callKey: CallKey => BypassedCallResults(
    //        Map(callKey -> (callKey.scope.outputs map { callOutput => callOutput.unqualifiedName -> JobOutput(WdlOptionalValue.none(callOutput.wdlType)) } toMap)))
    //      case declKey: DeclarationKey => BypassedDeclaration(declKey)
    //      case _ => throw new RuntimeException("Only calls and declarations might generate results when Bypassed")
    //    }
    throw new RuntimeException("Only calls and declarations might generate results when Bypassed")
  }

  def processRunnableDynamicDeclaration(declaration: DynamicDeclarationKey, data: WorkflowExecutionActorData) = {
    import lenthall.validation.ErrorOr._

    declaration.upstreamPorts.traverseValues(resolve(declaration, data)) map { lookup =>
      declaration.evaluate(lookup, data.expressionLanguageFunctions) match {
        case Valid(result) => self ! DeclarationEvaluationSucceededResponse(declaration, result)
        case Invalid(f) => self ! DeclarationEvaluationFailedResponse(declaration, new RuntimeException(f.toList.mkString(", ")))
      }
    } valueOr { f =>
      self ! DeclarationEvaluationFailedResponse(declaration, new RuntimeException(f.toList.mkString(", ")))
    }

    Success(WorkflowExecutionDiff(Map(declaration -> ExecutionStatus.Running)))
  }

  /**
    * Attempts to resolve an output port to a known value.
    * Curried for convenience.
    */
  private def resolve(jobKey: JobKey, data: WorkflowExecutionActorData)(outputPort: OutputPort): ErrorOr[WdlValue] = {
    data.outputStore.get(outputPort, jobKey.index).map(_.validNel) orElse {
      workflowDescriptor.defaultExpressions.get(outputPort) map {
        case expr if expr.inputs.isEmpty => expr.evaluateValue(Map.empty, data.expressionLanguageFunctions)
        case _ => "Cannot evaluate default expression with node dependencies".invalidNel
      }
    } getOrElse s"Can't find a value for ${outputPort.name}".invalidNel
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
    //    val scatterMap = conditionalKey.index flatMap { i =>
    //      // Will need update for nested scatters
    //      conditionalKey.scope.ancestry collectFirst { case s: Scatter => Map(s -> i) }
    //    } getOrElse Map.empty[Scatter, Int]
    //
    //    val lookup = conditionalKey.scope.lookupFunction(
    //      workflowDescriptor.knownValues,
    //      data.expressionLanguageFunctions,
    //      data.outputStore.fetchNodeOutputEntries,
    //      scatterMap
    //    )
    //
    //    conditionalKey.scope.condition.evaluate(lookup, data.expressionLanguageFunctions) map {
    //      case b: WdlBoolean =>
    //        val conditionalStatus = if (b.value) ExecutionStatus.Done else ExecutionStatus.Bypassed
    //        val result = WorkflowExecutionDiff(conditionalKey.populate(workflowDescriptor.knownValues) + (conditionalKey -> conditionalStatus))
    //        result
    //      case v: WdlValue => throw new RuntimeException(
    //        s"'if' condition must evaluate to a boolean but instead got ${v.wdlType.toWdlString}")
    //    }
    Failure(new Exception("BOOM"))
  }

  private def processRunnableScatter(scatterKey: ScatterKey, data: WorkflowExecutionActorData, bypassed: Boolean): Try[WorkflowExecutionDiff] = {
//    val lookup = scatterKey.scope.lookupFunction(
//      workflowDescriptor.knownValues,
//      data.expressionLanguageFunctions,
//      data.outputStore.fetchNodeOutputEntries
//    )
//
//    if (bypassed) {
//      Success(WorkflowExecutionDiff(scatterKey.populate(0, Map.empty) + (scatterKey -> ExecutionStatus.Bypassed)))
//    } else {
//      scatterKey.scope.collection.evaluate(lookup, data.expressionLanguageFunctions) map {
//        case WdlArrayLike(a) =>
//          WorkflowExecutionDiff(scatterKey.populate(a.value.size, workflowDescriptor.knownValues) + (scatterKey -> ExecutionStatus.Done))
//        case v: WdlValue => throw new RuntimeException(
//          s"Scatter collection must evaluate to an array but instead got ${v.wdlType.toWdlString}")
//      }
//    }
    Failure(new Exception("BOOM"))
  }

  private def processRunnableCollector(collector: CollectorKey, data: WorkflowExecutionActorData, isInBypassed: Boolean): Try[WorkflowExecutionDiff] = {

    val shards = data.executionStore.findCompletedShardsForOutput(collector)

    data.outputStore.generateCollectorOutput(collector, shards) match {
      case Failure(e) => Failure(new RuntimeException(s"Failed to collect output shards for call ${collector.tag}", e))
      case Success(outputs) =>
        val adjustedOutputs: CallOutputs = if (isInBypassed) {
          outputs map { output => (output._1, JobOutput(WdlOptionalValue.none(output._2.wdlValue.wdlType) )) }
        } else outputs
        self ! ScatterCollectionSucceededResponse(collector, adjustedOutputs)
        Success(WorkflowExecutionDiff(Map(collector -> ExecutionStatus.Starting)))
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

  case object RequestOutputStore extends WorkflowExecutionActorCommand

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

  private case class ScatterCollectionFailedResponse(collectorKey: CollectorKey, throwable: Throwable)

  private case class ScatterCollectionSucceededResponse(collectorKey: CollectorKey, outputs: CallOutputs)

  private case class DeclarationEvaluationSucceededResponse(declarationKey: ExpressionKey, value: WdlValue)

  private case object CheckRunnable

  private[execution] sealed trait BypassedScopeResults

  private case class BypassedCallResults(callOutputs: Map[CallKey, CallOutputs]) extends BypassedScopeResults
  private case class BypassedDeclaration(declaration: ExpressionKey) extends BypassedScopeResults

  private case class DeclarationEvaluationFailedResponse(declarationKey: ExpressionKey, reason: Throwable)

  case class SubWorkflowSucceededResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, outputs: CallOutputs)

  case class SubWorkflowFailedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, reason: Throwable)

  case class SubWorkflowAbortedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap)

  /**
    * Internal ADTs
    */
  case class ScatterKey(scatter: ScatterNode) extends JobKey {
    override val node = scatter
    override val index = None
    // When scatters are nested, this might become Some(_)
    override val attempt = 1
    override val tag = node.unqualifiedName

    /**
      * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
      *
      * @param count Number of ways to scatter the children.
      * @return ExecutionStore of scattered children.
      */
    def populate(count: Int, workflowCoercedInputs: WorkflowCoercedInputs): Map[JobKey, ExecutionStatus.Value] = {
      val keys = this.node.innerGraph.nodes flatMap {
        explode(_, count, workflowCoercedInputs)
      }
      keys map {
        _ -> ExecutionStatus.NotStarted
      } toMap
    }

    private def explode(scope: GraphNode, count: Int, workflowCoercedInputs: WorkflowCoercedInputs): Seq[JobKey] = {
      def makeCollectors(scope: GraphNode): Seq[CollectorKey] = scope match {
        case call: TaskCallNode => List(CollectorKey(call, scatter, count))
        case decl: ExpressionNode => List(CollectorKey(decl, scatter, count))
//        case i: If => i.children.flatMap(makeCollectors)
      }

      (scope match {
        case call: TaskCallNode => (0 until count) map { i => BackendJobDescriptorKey(call, Option(i), 1) }
//        case call: WdlWorkflowCall => (0 until count) map { i => SubWorkflowKey(call, Option(i), 1) }
//        case declaration: ExpressionNode => (0 until count) map { i => DeclarationKey(declaration, Option(i), workflowCoercedInputs) }
//        case conditional: If => (0 until count) map { i => ConditionalKey(conditional, Option(i)) }
//        case _: Scatter =>
//          throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
        case e =>
          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
      }) ++ makeCollectors(scope)
    }
  }

  // Represents a scatter collection for a call in the execution store
  case class CollectorKey(node: GraphNode, scatter: ScatterNode, scatterWidth: Int) extends JobKey {
    override val index = None
    override val attempt = 1
    override val tag = s"Collector-${node.unqualifiedName}"
  }

  case class SubWorkflowKey(node: WorkflowCallNode, index: ExecutionIndex, attempt: Int) extends CallKey {
    override val tag = s"SubWorkflow-${node.unqualifiedName}:${index.fromIndex}:$attempt"
  }

  case class ConditionalKey(ifScope: If, index: ExecutionIndex) extends JobKey {
    // TODO WOM: fixme
    override val node: GraphNode = null
    override val tag = node.unqualifiedName
    override val attempt = 1

    /**
      * Creates a sub-ExecutionStore with entries for each of the scoped children.
      *
      * @return ExecutionStore of scattered children.
      */
    def populate(workflowCoercedInputs: WorkflowCoercedInputs): Map[JobKey, ExecutionStatus.Value] = {
      //      scope.children map {
      //        keyify(_, workflowCoercedInputs) -> ExecutionStatus.NotStarted
      //      } toMap
      Map.empty
    }

//    /**
//      * Make a JobKey for all of the contained scopes.
//      */
//    private def keyify(scope: GraphNode, workflowCoercedInputs: WorkflowCoercedInputs): JobKey = {
//      scope match {
//        case call: TaskCallNode => BackendJobDescriptorKey(call, index, 1)
//        case call: WdlWorkflowCall => SubWorkflowKey(call, index, 1)
//        case declaration: Declaration => DeclarationKey(declaration, index, workflowCoercedInputs)
//        case i: If => ConditionalKey(i, index)
//        case scatter: Scatter if index.isEmpty => ScatterKey(scatter)
//        case _: Scatter =>
//          throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
//        case e =>
//          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported in an If block.")
//      }
//    }
  }

  object ExpressionKey {
    def apply(declaration: ExpressionNode, index: ExecutionIndex): ExpressionKey = {
      DynamicDeclarationKey(declaration, index)
    }
  }

  sealed trait ExpressionKey extends JobKey {
    override val node: ExpressionNode
    override val attempt = 1
    override val tag = s"Expression-${node.unqualifiedName}:${index.fromIndex}:$attempt"
  }

  case class DynamicDeclarationKey(node: ExpressionNode, index: ExecutionIndex) extends ExpressionKey {
    import lenthall.validation.ErrorOr._
    import lenthall.validation.Validation._

    lazy val inputs: Map[String, InputPort] = node.instantiatedExpression.inputMapping
    lazy val upstreamPorts: Map[String, OutputPort] = inputs map {
      case (key, input) => key -> input.upstream
    }

    def evaluate(lookup: Map[String, WdlValue], wdlFunctions: IoFunctionSet) = {
      node.instantiatedExpression.expression.evaluateValue(lookup, wdlFunctions) flatMap { node.womType.coerceRawValue(_).toErrorOr }
    }
  }

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
