package cromwell.engine.workflow.lifecycle.execution

import akka.actor._
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse, JobSucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.{AllBackendInitializationData, BackendJobDescriptorKey, JobExecutionMap}
import cromwell.core.Dispatcher._
import cromwell.core.ExecutionIndex._
import cromwell.core.ExecutionStatus._
import cromwell.core.WorkflowOptions.WorkflowFailureMode
import cromwell.core._
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{apply => _, _}
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, EngineLifecycleActorAbortedResponse}
import cromwell.engine.{ContinueWhilePossible, EngineWorkflowDescriptor}
import cromwell.services.metadata.MetadataService.{MetadataPutAcknowledgement, MetadataPutFailed}
import cromwell.util.StopAndLogSupervisor
import cromwell.webservice.EngineStatsActor
import lenthall.exception.ThrowableAggregation
import lenthall.util.TryUtil
import net.ceedubs.ficus.Ficus._
import wdl4s.values.{WdlArray, WdlBoolean, WdlOptionalValue, WdlValue, WdlString}
import org.apache.commons.lang3.StringUtils
import wdl4s.{Scope, _}

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class WorkflowExecutionActor(workflowDescriptor: EngineWorkflowDescriptor,
                                  serviceRegistryActor: ActorRef,
                                  jobStoreActor: ActorRef,
                                  subWorkflowStoreActor: ActorRef,
                                  callCacheReadActor: ActorRef,
                                  jobTokenDispenserActor: ActorRef,
                                  backendSingletonCollection: BackendSingletonCollection,
                                  initializationData: AllBackendInitializationData,
                                  restarting: Boolean)
  extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] with WorkflowLogging with CallMetadataHelper with StopAndLogSupervisor {
  
  implicit val ec = context.dispatcher
  
  override val workflowIdForLogging = workflowDescriptor.id
  override val workflowIdForCallMetadata = workflowDescriptor.id

  private val tag = s"WorkflowExecutionActor [UUID(${workflowDescriptor.id.shortString})]"
  private val MaxRetries = ConfigFactory.load().as[Option[Int]]("system.max-retries") match {
    case Some(value) => value
    case None =>
      workflowLogger.warn(s"Failed to load the max-retries value from the configuration. Defaulting back to a value of '$DefaultMaxRetriesFallbackValue'.")
      DefaultMaxRetriesFallbackValue
  }
  
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
      outputStore = OutputStore.empty
    )
  )

  when(WorkflowExecutionPendingState) {
    case Event(ExecuteWorkflowCommand, stateData) =>
      val data = startRunnableScopes(stateData)
      goto(WorkflowExecutionInProgressState) using data
  }

  when(WorkflowExecutionInProgressState) {
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
    case Event(JobSucceededResponse(jobKey, returnCode, callOutputs, _, _), stateData) =>
      pushSuccessfulCallMetadata(jobKey, returnCode, callOutputs)
      handleCallSuccessful(jobKey, callOutputs, stateData, Map.empty)
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
      handleDeclarationEvaluationSuccessful(declKey, WdlOptionalValue.none(declKey.scope.wdlType), stateData)

      // Failure
        // Initialization
    case Event(JobInitializationFailed(jobKey, reason), stateData) =>
      pushFailedCallMetadata(jobKey, None, reason, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, reason, Map.empty)
      // Job Non Retryable
    case Event(JobFailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedCallMetadata(jobKey, returnCode, reason, retryableFailure = false)
      handleNonRetryableFailure(stateData, jobKey, reason, Map.empty)
      // Job Retryable
    case Event(JobFailedRetryableResponse(jobKey, reason, returnCode), stateData) =>
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
      handleCallAborted(stateData, jobKey, Map.empty)
    case Event(SubWorkflowAbortedResponse(jobKey, executedKeys), stateData) =>
      handleCallAborted(stateData, jobKey, executedKeys)
    case Event(SubWorkflowSucceededResponse(subKey, executedKeys, _), stateData) =>
      handleCallAborted(stateData, subKey, executedKeys)
    case Event(JobSucceededResponse(jobKey, returnCode, callOutputs, _, _), stateData) =>
      handleCallAborted(stateData, jobKey, Map.empty)
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
      pushFailedCallMetadata(jobKey, None, reason, retryableFailure = false)
      stay
    case Event(JobFailedNonRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedCallMetadata(jobKey, returnCode, reason, retryableFailure = false)
      stay
    case Event(JobFailedRetryableResponse(jobKey, reason, returnCode), stateData) =>
      pushFailedCallMetadata(jobKey, returnCode, reason, retryableFailure = true)
      stay
    case Event(JobSucceededResponse(jobKey, returnCode, callOutputs, _, _), stateData) =>
      pushSuccessfulCallMetadata(jobKey, returnCode, callOutputs)
      stay
  }


  def handleTerminated(actorRef: ActorRef) = {
    // Both of these Should Never Happen (tm), assuming the state data is set correctly on EJEA creation.
    // If they do, it's a big programmer error and the workflow execution fails.
    val jobKey = stateData.engineCallExecutionActors.getOrElse(actorRef, throw new RuntimeException("Programmer Error: An EJEA has terminated but was not assigned a jobKey"))
    val jobStatus = stateData.executionStore.store.getOrElse(jobKey, throw new RuntimeException("Programmer Error: An EJEA representing a jobKey which this workflow is not running has sent up a terminated message."))

    if (!jobStatus.isTerminal) {
      val terminationException = getFailureCause(actorRef) match {
        case Some(e) => new RuntimeException("Unexpected failure in EJEA.", e)
        case None => new RuntimeException("Unexpected failure in EJEA (root cause not captured).")
      }
      self ! JobFailedNonRetryableResponse(jobKey, terminationException, None)
    }

    stay
  }

  whenUnhandled {
    case Event(Terminated(actorRef), stateData) => handleTerminated(actorRef) using stateData.removeEngineJobExecutionActor(actorRef)
    case Event(MetadataPutFailed(action, error), _) =>
      // Do something useful here??
      workflowLogger.warn(s"$tag Put failed for Metadata action $action", error)
      stay()
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    case Event(EngineLifecycleActorAbortCommand, stateData) =>
      if (stateData.hasRunningActors) {
        log.info(s"$tag: Abort received. " +
          s"Aborting ${stateData.backendJobExecutionActors.size} Job Execution Actors" +
          s"and ${stateData.subWorkflowExecutionActors.size} Sub Workflow Execution Actors"
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
  
  private def handleDeclarationEvaluationFailure(declarationKey: DeclarationKey, reason: Throwable, stateData: WorkflowExecutionActorData) = {
    handleExecutionFailure(declarationKey, stateData, reason, Map.empty)
  }
  
  private def handleExecutionFailure(failedJobKey: JobKey, data: WorkflowExecutionActorData, reason: Throwable, jobExecutionMap: JobExecutionMap) = {
    val newData = data.executionFailed(failedJobKey)
    
    if (workflowDescriptor.getWorkflowOption(WorkflowFailureMode).contains(ContinueWhilePossible.toString)) {
      newData.workflowCompletionStatus match {
        case Some(completionStatus) if completionStatus == Failed =>
          context.parent ! WorkflowExecutionFailedResponse(newData.jobExecutionMap, reason)
          goto(WorkflowExecutionFailedState) using newData
        case _ =>
          stay() using startRunnableScopes(newData)
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

     val (response, finalState) = workflowDescriptor.workflow.evaluateOutputs(
      workflowDescriptor.knownValues,
      data.expressionLanguageFunctions,
      data.outputStore.fetchNodeOutputEntries
    ) map { workflowOutputs =>
       workflowLogger.info(
         s"""Workflow ${workflowDescriptor.workflow.unqualifiedName} complete. Final Outputs:
             |${workflowOutputs.stripLarge.toJson.prettyPrint}""".stripMargin
       )
       pushWorkflowOutputMetadata(workflowOutputs)
       (WorkflowExecutionSucceededResponse(data.jobExecutionMap, workflowOutputs mapValues JobOutput.apply), WorkflowExecutionSuccessfulState)
    } recover {
       case ex =>
         (WorkflowExecutionFailedResponse(data.jobExecutionMap, ex), WorkflowExecutionFailedState)
    } get
    
    context.parent ! response
    goto(finalState) using data
  }

  private def handleRetryableFailure(jobKey: BackendJobDescriptorKey, reason: Throwable, returnCode: Option[Int]) = {
    // We start with index 1 for #attempts, hence invariant breaks only if jobKey.attempt > MaxRetries
    if (jobKey.attempt <= MaxRetries) {
      val newJobKey = jobKey.copy(attempt = jobKey.attempt + 1)
      workflowLogger.info(s"Retrying job execution for ${newJobKey.tag}")
      /*  Currently, we update the status of the old key to Preempted, and add a new entry (with the #attempts incremented by 1)
        * to the execution store with status as NotStarted. This allows startRunnableCalls to re-execute this job */
      val executionDiff = WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Preempted, newJobKey -> ExecutionStatus.NotStarted))
      val newData = stateData.mergeExecutionDiff(executionDiff).removeCallExecutionActor(jobKey)
      stay() using startRunnableScopes(newData)
    } else {
      workflowLogger.warn(s"Exhausted maximum number of retries for job ${jobKey.tag}. Failing.")
      goto(WorkflowExecutionFailedState) using stateData.mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed))).removeCallExecutionActor(jobKey)
    }
  }

  private def handleCallSuccessful(jobKey: JobKey, outputs: CallOutputs, data: WorkflowExecutionActorData, jobExecutionMap: JobExecutionMap) = {
    handleExecutionSuccess(data.callExecutionSuccess(jobKey, outputs).addExecutions(jobExecutionMap))
  }
  
  private def handleDeclarationEvaluationSuccessful(key: DeclarationKey, value: WdlValue, data: WorkflowExecutionActorData) = {
    handleExecutionSuccess(data.declarationEvaluationSuccess(key, value))
  }

  private def handleCallBypassed(callOutputs: Map[CallKey, CallOutputs], data: WorkflowExecutionActorData) = {
    def foldFunc(d: WorkflowExecutionActorData, output: (CallKey, CallOutputs)) = d.callExecutionSuccess(output._1, output._2)

    val updatedData = callOutputs.foldLeft(data)(foldFunc)
    handleExecutionSuccess(updatedData)
  }

  private def handleExecutionSuccess(data: WorkflowExecutionActorData) = {
    data.workflowCompletionStatus match {
      case Some(ExecutionStatus.Done) =>
        handleWorkflowSuccessful(data)
      case Some(sts) =>
        context.parent ! WorkflowExecutionFailedResponse(data.jobExecutionMap, new Exception("One or more jobs failed in fail-slow mode"))
        goto(WorkflowExecutionFailedState) using data
      case _ =>
        stay() using startRunnableScopes(data)
    }
  }
  
  private def handleCallAborted(data: WorkflowExecutionActorData, jobKey: JobKey, jobExecutionMap: JobExecutionMap) = {
    workflowLogger.info(s"$tag job aborted: ${jobKey.tag}")
    val newStateData = data.removeCallExecutionActor(jobKey).addExecutions(jobExecutionMap)
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
  @tailrec
  private def startRunnableScopes(data: WorkflowExecutionActorData): WorkflowExecutionActorData = {
    val runnableScopes = data.executionStore.runnableScopes
    val runnableCalls = runnableScopes.view collect { case k if k.scope.isInstanceOf[Call] => k } sortBy { k =>
      (k.scope.fullyQualifiedName, k.index.getOrElse(-1)) } map { _.tag }
    if (runnableCalls.nonEmpty) workflowLogger.info("Starting calls: " + runnableCalls.mkString(", "))

    // Each process returns a Try[WorkflowExecutionDiff], which, upon success, contains potential changes to be made to the execution store.
    val executionDiffs = runnableScopes map {
      case k: CallKey if isInBypassedScope(k, data) => processBypassedScope(k, data)
      case k: DeclarationKey if isInBypassedScope(k, data) => processBypassedScope(k, data)
      case k: BackendJobDescriptorKey => processRunnableJob(k, data)
      case k: ScatterKey => processRunnableScatter(k, data, isInBypassedScope(k, data))
      case k: ConditionalKey => processRunnableConditional(k, data)
      case k: CollectorKey => processRunnableCollector(k, data, isInBypassedScope(k, data))
      case k: SubWorkflowKey => processRunnableSubWorkflow(k, data)
      case k: StaticDeclarationKey => processRunnableStaticDeclaration(k)
      case k: DynamicDeclarationKey => processRunnableDynamicDeclaration(k, data)
      case k =>
        val exception = new UnsupportedOperationException(s"Unknown entry in execution store: ${k.tag}")
        self ! JobInitializationFailed(k, exception)
        Failure(exception)
    }

    TryUtil.sequence(executionDiffs) match {
      case Success(diffs) =>
        // Update the metadata for the jobs we just sent to EJEAs (they'll start off queued up waiting for tokens):
        pushQueuedCallMetadata(diffs)
        if (diffs.exists(_.containsNewEntry)) {
          val newData = data.mergeExecutionDiffs(diffs)
          startRunnableScopes(newData)
        } else {
          val result = data.mergeExecutionDiffs(diffs)
          result
        }
      case Failure(e) => throw new RuntimeException("Unexpected engine failure", e)
    }
  }

  private def isInBypassedScope(jobKey: JobKey, data: WorkflowExecutionActorData) = {
    // This is hugely ripe for optimization, if it becomes a bottleneck
    def isBypassedConditional(conditional: If, executionStore: ExecutionStore): Boolean = {
      executionStore.store.exists { case (jobKeyUnderExamination, executionStatus) =>
        if (jobKeyUnderExamination.isInstanceOf[ConditionalKey]) {
          if (jobKeyUnderExamination.scope.fullyQualifiedName.equals(conditional.fullyQualifiedName) && jobKeyUnderExamination.index.equals(jobKey.index)) {
            executionStatus.equals(ExecutionStatus.Bypassed)
          } else false
        } else false
      }
    }

    val result = jobKey.scope.ancestry.exists {
      case i: If => isBypassedConditional(i, data.executionStore)
      case _ => false
    }
    result
  }

  def processBypassedScope(jobKey: JobKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    self ! bypassedScopeResults(jobKey)
    Success(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Running)))
  }

  def bypassedScopeResults(jobKey: JobKey): BypassedScopeResults = jobKey match {
    case callKey: CallKey => BypassedCallResults(
      Map(callKey -> (callKey.scope.outputs map { callOutput => callOutput.unqualifiedName -> JobOutput(WdlOptionalValue.none(callOutput.wdlType)) } toMap)))
    case declKey: DeclarationKey => BypassedDeclaration(declKey)
    case _ => throw new RuntimeException("Only calls and declarations might generate results when Bypassed")
  }

  def processRunnableStaticDeclaration(declaration: StaticDeclarationKey) = {
    self ! DeclarationEvaluationSucceededResponse(declaration, declaration.value)
    Success(WorkflowExecutionDiff(Map(declaration -> ExecutionStatus.Running)))
  }
  
  def processRunnableDynamicDeclaration(declaration: DynamicDeclarationKey, data: WorkflowExecutionActorData) = {
    val scatterMap = declaration.index flatMap { i =>
      // Will need update for nested scatters
      declaration.scope.ancestry collectFirst { case s: Scatter => Map(s -> i) }
    } getOrElse Map.empty[Scatter, Int]

    val lookup = declaration.scope.lookupFunction(
      workflowDescriptor.knownValues,
      data.expressionLanguageFunctions,
      data.outputStore.fetchNodeOutputEntries,
      scatterMap
    )

    declaration.requiredExpression.evaluate(lookup, data.expressionLanguageFunctions) match {
      case Success(result) => self ! DeclarationEvaluationSucceededResponse(declaration, result)
      case Failure(ex) => self ! DeclarationEvaluationFailedResponse(declaration, ex)
    }

    Success(WorkflowExecutionDiff(Map(declaration -> ExecutionStatus.Running)))
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
              self, jobKey, data, factory, initializationData.get(backendName), restarting, serviceRegistryActor,
              jobStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingleton, backendName, workflowDescriptor.callCachingMode)
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
      SubWorkflowExecutionActor.props(key, data, backendFactories, serviceRegistryActor, jobStoreActor, subWorkflowStoreActor,
        callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection, initializationData, restarting),
      s"SubWorkflowExecutionActor-${key.tag}"
    )

    context watch sweaRef
    pushNewCallMetadata(key, None)
    sweaRef ! SubWorkflowExecutionActor.Execute
    
    Success(WorkflowExecutionDiff(executionStoreChanges = Map(key -> ExecutionStatus.QueuedInCromwell),
      engineJobExecutionActorAdditions = Map(sweaRef -> key)))
  }

  private def processRunnableConditional(conditionalKey: ConditionalKey, data: WorkflowExecutionActorData): Try[WorkflowExecutionDiff] = {
    val scatterMap = conditionalKey.index flatMap { i =>
      // Will need update for nested scatters
      conditionalKey.scope.ancestry collectFirst { case s: Scatter => Map(s -> i) }
    } getOrElse Map.empty[Scatter, Int]

    val lookup = conditionalKey.scope.lookupFunction(
      workflowDescriptor.knownValues,
      data.expressionLanguageFunctions,
      data.outputStore.fetchNodeOutputEntries,
      scatterMap
    )

    conditionalKey.scope.condition.evaluate(lookup, data.expressionLanguageFunctions) map {
      case b: WdlBoolean =>
        val conditionalStatus = if (b.value) ExecutionStatus.Done else ExecutionStatus.Bypassed
        val result = WorkflowExecutionDiff(conditionalKey.populate(workflowDescriptor.knownValues) + (conditionalKey -> conditionalStatus))
        result
      case v: WdlValue => throw new RuntimeException("'if' condition must evaluate to a boolean")
    }
  }

  private def processRunnableScatter(scatterKey: ScatterKey, data: WorkflowExecutionActorData, bypassed: Boolean): Try[WorkflowExecutionDiff] = {
    val lookup = scatterKey.scope.lookupFunction(
      workflowDescriptor.knownValues,
      data.expressionLanguageFunctions,
      data.outputStore.fetchNodeOutputEntries
    )

    if (bypassed) {
      Success(WorkflowExecutionDiff(scatterKey.populate(0, Map.empty) + (scatterKey -> ExecutionStatus.Bypassed)))
    } else {
      scatterKey.scope.collection.evaluate(lookup, data.expressionLanguageFunctions) map {
        case a: WdlArray =>
          WorkflowExecutionDiff(scatterKey.populate(a.value.size, workflowDescriptor.knownValues) + (scatterKey -> ExecutionStatus.Done))
        case v: WdlValue => throw new RuntimeException("Scatter collection must evaluate to an array")
      }
    }
  }

  private def processRunnableCollector(collector: CollectorKey, data: WorkflowExecutionActorData, isInBypassed: Boolean): Try[WorkflowExecutionDiff] = {
    def shards(collectorKey: CollectorKey) = data.executionStore.findShardEntries(collectorKey) collect {
      case (k: CallKey, v) if v.isDoneOrBypassed => k
      case (k: DynamicDeclarationKey, v) if v.isDoneOrBypassed => k
    }
    data.outputStore.generateCollectorOutput(collector, shards(collector)) match {
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
  
  private case class DeclarationEvaluationSucceededResponse(declarationKey: DeclarationKey, value: WdlValue)

  private[execution] sealed trait BypassedScopeResults

  private case class BypassedCallResults(callOutputs: Map[CallKey, CallOutputs]) extends BypassedScopeResults
  private case class BypassedDeclaration(declaration: DeclarationKey) extends BypassedScopeResults

  private case class DeclarationEvaluationFailedResponse(declarationKey: DeclarationKey, reason: Throwable)

  case class SubWorkflowSucceededResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, outputs: CallOutputs)

  case class SubWorkflowFailedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap, reason: Throwable)

  case class SubWorkflowAbortedResponse(key: SubWorkflowKey, jobExecutionMap: JobExecutionMap)

  /**
    * Internal ADTs
    */
  case class ScatterKey(scatter: Scatter) extends JobKey {
    override val scope = scatter
    override val index = None
    // When scatters are nested, this might become Some(_)
    override val attempt = 1
    override val tag = scope.unqualifiedName

    /**
      * Creates a sub-ExecutionStore with Starting entries for each of the scoped children.
      *
      * @param count Number of ways to scatter the children.
      * @return ExecutionStore of scattered children.
      */
    def populate(count: Int, workflowCoercedInputs: WorkflowCoercedInputs): Map[JobKey, ExecutionStatus.Value] = {
      val keys = this.scope.children flatMap {
        explode(_, count, workflowCoercedInputs)
      }
      keys map {
        _ -> ExecutionStatus.NotStarted
      } toMap
    }

    private def explode(scope: Scope, count: Int, workflowCoercedInputs: WorkflowCoercedInputs): Seq[JobKey] = {
      def makeCollectors(scope: Scope): Seq[CollectorKey] = scope match {
        case call: Call => List(CollectorKey(call, scatter, count))
        case decl: Declaration => List(CollectorKey(decl, scatter, count))
        case i: If => i.children.flatMap(makeCollectors(_))
      }

      (scope match {
        case call: TaskCall => (0 until count) map { i => BackendJobDescriptorKey(call, Option(i), 1) }
        case call: WorkflowCall => (0 until count) map { i => SubWorkflowKey(call, Option(i), 1) }
        case declaration: Declaration => (0 until count) map { i => DeclarationKey(declaration, Option(i), workflowCoercedInputs) }
        case conditional: If => (0 until count) map { i => ConditionalKey(conditional, Option(i)) }
        case scatter: Scatter =>
          throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
        case e =>
          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported.")
      }) ++ makeCollectors(scope)
    }
  }

  // Represents a scatter collection for a call in the execution store
  case class CollectorKey(scope: Scope with GraphNode, scatter: Scatter, scatterWidth: Int) extends JobKey {
    override val index = None
    override val attempt = 1
    override val tag = s"Collector-${scope.unqualifiedName}"
  }

  case class SubWorkflowKey(scope: WorkflowCall, index: ExecutionIndex, attempt: Int) extends CallKey {
    override val tag = s"SubWorkflow-${scope.unqualifiedName}:${index.fromIndex}:$attempt"
  }

  case class ConditionalKey(scope: If, index: ExecutionIndex) extends JobKey {

    override val tag = scope.unqualifiedName
    override val attempt = 1

    /**
      * Creates a sub-ExecutionStore with entries for each of the scoped children.
      *
      * @return ExecutionStore of scattered children.
      */
    def populate(workflowCoercedInputs: WorkflowCoercedInputs): Map[JobKey, ExecutionStatus.Value] = {
      scope.children map {
        keyify(_, workflowCoercedInputs) -> ExecutionStatus.NotStarted
      } toMap
    }

    /**
      * Make a JobKey for all of the contained scopes.
      */
    private def keyify(scope: Scope, workflowCoercedInputs: WorkflowCoercedInputs): JobKey = {
      scope match {
        case call: TaskCall => BackendJobDescriptorKey(call, index, 1)
        case call: WorkflowCall => SubWorkflowKey(call, index, 1)
        case declaration: Declaration => DeclarationKey(declaration, index, workflowCoercedInputs)
        case i: If => ConditionalKey(i, index)
        case scatter: Scatter if index.isEmpty => ScatterKey(scatter)
        case _: Scatter =>
          throw new UnsupportedOperationException("Nested Scatters are not supported (yet) ... but you might try a sub workflow to achieve the same effect!")
        case e =>
          throw new UnsupportedOperationException(s"Scope ${e.getClass.getName} is not supported in an If block.")
      }
    }
  }

  object DeclarationKey {
    def apply(declaration: Declaration, index: ExecutionIndex, inputs: WorkflowCoercedInputs): DeclarationKey = {
      inputs.find(_._1 == declaration.fullyQualifiedName) match {
        case Some((_, value)) => StaticDeclarationKey(declaration, index, value)
        case None => declaration.expression map { expression =>
          DynamicDeclarationKey(declaration, index, expression)
        } getOrElse {
          throw new RuntimeException(s"Found a declaration ${declaration.fullyQualifiedName} without expression and without input value. This should have been a validation error.")
        }
      }
    }
  }
    
  sealed trait DeclarationKey extends JobKey {
    override val scope: Declaration
    override val attempt = 1
    override val tag = s"Declaration-${scope.unqualifiedName}:${index.fromIndex}:$attempt"
  }
  
  case class StaticDeclarationKey(scope: Declaration, index: ExecutionIndex, value: WdlValue) extends DeclarationKey
  
  case class DynamicDeclarationKey(scope: Declaration, index: ExecutionIndex, requiredExpression: WdlExpression) extends DeclarationKey

  case class WorkflowExecutionException[T <: Throwable](exceptions: NonEmptyList[T]) extends ThrowableAggregation {
    override val throwables = exceptions.toList
    override val exceptionContext = s"WorkflowExecutionActor"
  }

  private lazy val DefaultMaxRetriesFallbackValue = 10

  def props(workflowDescriptor: EngineWorkflowDescriptor,
            serviceRegistryActor: ActorRef,
            jobStoreActor: ActorRef,
            subWorkflowStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection,
            initializationData: AllBackendInitializationData,
            restarting: Boolean): Props = {
    Props(WorkflowExecutionActor(workflowDescriptor, serviceRegistryActor, jobStoreActor, subWorkflowStoreActor,
      callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection, initializationData, restarting)).withDispatcher(EngineDispatcher)
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
