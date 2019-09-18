package cromwell.engine.workflow.lifecycle.execution.job

import akka.actor.SupervisorStrategy.{Escalate, Stop}
import akka.actor.{ActorInitializationException, ActorRef, LoggingFSM, OneForOneStrategy, Props}
import cromwell.backend.BackendCacheHitCopyingActor.{CopyOutputsCommand, CopyingOutputsFailedResponse}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend._
import cromwell.backend.standard.StandardInitializationData
import cromwell.backend.standard.callcaching.BlacklistCache
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.CallCachingEntry
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.instrumentation.JobInstrumentation
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestValueStore
import cromwell.engine.workflow.lifecycle.execution._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed}
import cromwell.engine.workflow.lifecycle.execution.job.preparation.{CallPreparation, JobPreparationActor}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.engine.workflow.lifecycle.{EngineLifecycleActorAbortCommand, TimedFSM}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore._
import cromwell.services.EngineServicesStore
import cromwell.services.metadata.CallMetadataKeys.CallCachingKeys
import cromwell.services.metadata.{CallMetadataKeys, MetadataKey}
import cromwell.webservice.EngineStatsActor

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

class EngineJobExecutionActor(replyTo: ActorRef,
                              jobDescriptorKey: BackendJobDescriptorKey,
                              workflowDescriptor: EngineWorkflowDescriptor,
                              backendLifecycleActorFactory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              val serviceRegistryActor: ActorRef,
                              ioActor: ActorRef,
                              jobStoreActor: ActorRef,
                              callCacheReadActor: ActorRef,
                              callCacheWriteActor: ActorRef,
                              workflowDockerLookupActor: ActorRef,
                              jobTokenDispenserActor: ActorRef,
                              backendSingletonActor: Option[ActorRef],
                              callCachingMode: CallCachingMode,
                              command: BackendJobExecutionActorCommand,
                              fileHashCachingActor: Option[ActorRef],
                              blacklistCache: Option[BlacklistCache]) extends LoggingFSM[EngineJobExecutionActorState, EJEAData]
  with WorkflowLogging with CallMetadataHelper with JobInstrumentation with TimedFSM[EngineJobExecutionActorState] {

  override val workflowIdForLogging = workflowDescriptor.possiblyNotRootWorkflowId
  override val rootWorkflowIdForLogging = workflowDescriptor.rootWorkflowId
  override val workflowIdForCallMetadata = workflowDescriptor.id

  val jobStartTime = System.currentTimeMillis()

  override val supervisorStrategy = OneForOneStrategy() {
    // If an actor fails to initialize, send the exception to self before stopping it so we can fail the job properly
    case e: ActorInitializationException =>
      self ! e
      Stop
    case t =>
      super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }

  val jobTag = s"${workflowIdForLogging.shortString}:${jobDescriptorKey.call.fullyQualifiedName}:${jobDescriptorKey.index.fromIndex}:${jobDescriptorKey.attempt}"
  val tag = s"EJEA_$jobTag"

  //noinspection ActorMutableStateInspection
  // There's no need to check for a cache hit again if we got preempted, or if there's no result copying actor defined
  // NB: this can also change (e.g. if we have a HashError we just force this to CallCachingOff)
  private[execution] var effectiveCallCachingMode = {
    if (backendLifecycleActorFactory.fileHashingActorProps.isEmpty) CallCachingOff
    else if (backendLifecycleActorFactory.cacheHitCopyingActorProps.isEmpty || jobDescriptorKey.attempt > 1) {
      callCachingMode.withoutRead
    } else callCachingMode
  }

  // For tests:
  private[execution] def checkEffectiveCallCachingMode = effectiveCallCachingMode

  private val effectiveCallCachingKey = CallCachingKeys.EffectiveModeKey
  private val callCachingReadResultMetadataKey = CallCachingKeys.ReadResultMetadataKey
  private val callCachingHitResultMetadataKey = CallCachingKeys.HitResultMetadataKey
  private val callCachingAllowReuseMetadataKey = CallCachingKeys.AllowReuseMetadataKey
  private val callCachingHashes = CallCachingKeys.HashesKey

  private val callCachingOptionsOption = effectiveCallCachingMode match {
    case callCachingActivity: CallCachingActivity => Option(callCachingActivity.options)
    case _ => None
  }

  private val invalidationRequired = callCachingOptionsOption.exists(_.invalidateBadCacheResults)

  private val callCachePathPrefixes = for {
    callCachingOptions <- callCachingOptionsOption
    workflowOptionPrefixes <- callCachingOptions.workflowOptionCallCachePrefixes
    d <- initializationData collect { case d: StandardInitializationData => d }
    rootPrefix = d.workflowPaths.callCacheRootPrefix
  } yield CallCachePathPrefixes(rootPrefix, workflowOptionPrefixes.toList)

  implicit val ec: ExecutionContext = context.dispatcher

  override def preStart() = {
    log.debug(s"$tag: $effectiveCallCachingKey: $effectiveCallCachingMode")
  }

  startWith(Pending, NoData)
  //noinspection ActorMutableStateInspection
  private var eventList: Seq[ExecutionEvent] = Seq(ExecutionEvent(stateName.toString))

  override def onTimedTransition(from: EngineJobExecutionActorState, to: EngineJobExecutionActorState, duration: FiniteDuration) = {
    // Send to StatsD
    recordExecutionStepTiming(from.toString, duration)
  }

  // When Pending, the FSM always has NoData
  when(Pending) {
    case Event(Execute, NoData) =>
      requestExecutionToken()
      goto(RequestingExecutionToken)
  }

  when(RequestingExecutionToken) {
    case Event(JobExecutionTokenDispensed, NoData) =>
      replyTo ! JobStarting(jobDescriptorKey)
      if (restarting) {
        val jobStoreKey = jobDescriptorKey.toJobStoreKey(workflowIdForLogging)
        jobStoreActor ! QueryJobCompletion(jobStoreKey, jobDescriptorKey.call.outputPorts.toSeq)
        goto(CheckingJobStore)
      } else {
        requestValueStore()
      }
  }

  // When CheckingJobStore, the FSM always has NoData
  when(CheckingJobStore) {
    case Event(JobNotComplete, NoData) =>
      checkCacheEntryExistence()
    case Event(JobComplete(jobResult), NoData) =>
      respondAndStop(jobResult.toBackendJobResponse(jobDescriptorKey))
    case Event(f: JobStoreReadFailure, NoData) =>
      writeCallCachingModeToMetadata()
      log.error(f.reason, "{}: Error reading from JobStore", tag)
      // Escalate
      throw new RuntimeException(f.reason)
  }

  // When we're restarting but the job store says the job is not complete.
  // This is to cover for the case where Cromwell was stopped after writing Cache Info to the DB but before
  // writing to the JobStore. In that case, we can re-use the cached results and save them to the job store directly.
  // Note that this checks that *this* particular job has a cache entry, not that there is *a* cache hit (possibly from another job)
  // There's no need to copy the outputs because they're already *this* job's outputs
  when(CheckingCacheEntryExistence) {
    // There was already a cache entry for this job
    case Event(join: CallCachingJoin, NoData) =>
      Try(join.toJobSuccess(jobDescriptorKey, backendLifecycleActorFactory.pathBuilders(initializationData))).map({ jobSuccess =>
        // We can't create a CallCacheHashes to give to the SucceededResponseData here because it involves knowledge of
        // which hashes are file hashes and which are not. We can't know that (nor do we care) when pulling them from the
        // database. So instead manually publish the hashes here.
        publishHashResultsToMetadata(Option(Success(join.callCacheHashes)))
        saveJobCompletionToJobStore(SucceededResponseData(jobSuccess, None))
      }).recover({
        case f =>
          // If for some reason the above fails, fail the job cleanly
          saveJobCompletionToJobStore(FailedResponseData(JobFailedNonRetryableResponse(jobDescriptorKey, f, None), None))
      }).get
    // No cache entry for this job - keep going
    case Event(NoCallCacheEntry(_), NoData) =>
      requestValueStore()
    case Event(CacheResultLookupFailure(reason), NoData) =>
      log.error(reason, "{}: Failure checking for cache entry existence: {}. Attempting to resume job anyway.", jobTag, reason.getMessage)
      requestValueStore()
  }

  /*
  * ! Hot Potato Warning !
  * We ask explicitly for the output store so we can use it on the fly and more importantly not store it as a
  * variable in this actor, which would prevent it from being garbage collected for the duration of the
  * job and would lead to memory leaks.
  */
  when(WaitingForValueStore) {
    case Event(valueStore: ValueStore, NoData) => prepareJob(valueStore)
  }

  // When PreparingJob, the FSM always has NoData
  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, bjeaProps), NoData) =>
      val updatedData = ResponsePendingData(jobDescriptor, bjeaProps)
      effectiveCallCachingMode match {
        case activity: CallCachingActivity if activity.readFromCache => handleReadFromCacheOn(jobDescriptor, activity, updatedData)
        case activity: CallCachingActivity => handleReadFromCacheOff(jobDescriptor, activity, updatedData)
        case CallCachingOff => runJob(updatedData)
      }
    case Event(CallPreparationFailed(jobKey: BackendJobDescriptorKey, throwable), NoData) =>
      respondAndStop(JobFailedNonRetryableResponse(jobKey, throwable, None))
  }

  when(CheckingCallCache) {
    case Event(CacheMiss, data: ResponsePendingData) =>
      writeToMetadata(Map(
        callCachingHitResultMetadataKey -> false,
        callCachingReadResultMetadataKey -> "Cache Miss"))
      log.debug("Cache miss for job {}", jobTag)
      runJob(data)
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
      addHashesAndStay(data, hashes)
    case Event(hit: CacheHit, data: ResponsePendingData) =>
      fetchCachedResults(hit.cacheResultId, data.withCacheHit(hit))
    case Event(HashError(t), data: ResponsePendingData) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Hashing Error: ${t.getMessage}"))
      disableCallCaching(Option(t))
      runJob(data)
  }

  when(FetchingCachedOutputsFromDatabase) {
    case Event(
    CachedOutputLookupSucceeded(womValueSimpletons, jobDetritus, returnCode, cacheResultId, cacheHitDetails),
    data@ResponsePendingData(_, _, _, _, Some(ejeaCacheHit), _, _),
    ) =>
      if (cacheResultId != ejeaCacheHit.hit.cacheResultId) {
        // Sanity check: was this the right set of results (a false here is a BAD thing!):
        log.error(s"Received incorrect call cache results from FetchCachedResultsActor. Expected ${ejeaCacheHit.hit} but got $cacheResultId. Running job")
        // Treat this like the "CachedOutputLookupFailed" event:
        runJob(data)
      } else {
        writeToMetadata(Map(
          callCachingHitResultMetadataKey -> true,
          callCachingReadResultMetadataKey -> s"Cache Hit: $cacheHitDetails"))
        log.debug("Cache hit for {}! Fetching cached result {}", jobTag, cacheResultId)
        makeBackendCopyCacheHit(womValueSimpletons, jobDetritus, returnCode, data, cacheResultId, ejeaCacheHit.hitNumber) using data.withCacheDetails(cacheHitDetails)
      }
    case Event(CachedOutputLookupFailed(_, error), data: ResponsePendingData) =>
      log.warning("Can't fetch a list of cached outputs to copy for {} due to {}. Running job.", jobTag, error)
      runJob(data)
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
      addHashesAndStay(data, hashes)
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCacheWrite(t)
      // Can't write hashes for this job, but continue to wait for the lookup response.
      stay using data.copy(hashes = Option(Failure(t)))
  }

  when(BackendIsCopyingCachedOutputs) {
    // Backend copying response:
    case Event(
    response: JobSucceededResponse,
    data@ResponsePendingData(_, _, Some(Success(hashes)), _, _, _, _),
    ) =>
      logCacheHitSuccess(data)
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: JobSucceededResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      logCacheHitSuccess(data)
      // Wait for the CallCacheHashes
      stay using data.withSuccessResponse(response)
    case Event(response: JobSucceededResponse, data: ResponsePendingData) => // bad hashes or cache write off
      logCacheHitSuccess(data)
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
    case Event(
    CopyingOutputsFailedResponse(_, cacheCopyAttempt, throwable),
    data@ResponsePendingData(_, _, _, _, Some(cacheHit), _, _)
    ) if cacheCopyAttempt == cacheHit.hitNumber =>
      invalidateCacheHitAndTransition(cacheHit, data, throwable)

    // Hashes arrive:
    case Event(hashes: CallCacheHashes, data: SucceededResponseData) =>
      saveCacheResults(hashes, data)
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
      addHashesAndStay(data, hashes)

    // Hash error occurs:
    case Event(HashError(t), data: SucceededResponseData) =>
      disableCacheWrite(t)
      saveJobCompletionToJobStore(data.copy(hashes = Option(Failure(t))))
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCacheWrite(t)
      // Can't write hashes for this job, but continue to wait for the copy response.
      stay using data.copy(hashes = Option(Failure(t)))
  }

  when(InvalidatingCacheEntry) {
    case Event(response: CallCacheInvalidatedResponse, data: ResponsePendingData) =>
      handleCacheInvalidatedResponse(response, data)

    // Hashes arrive:
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
      addHashesAndStay(data, hashes)

    // Hash error occurs:
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCacheWrite(t)
      // Can't write hashes for this job, but continue to wait for the copy response.
      stay using data.copy(hashes = Option(Failure(t)))
  }

  // Handles JobSucceededResponse messages
  val jobSuccessHandler: StateFunction = {
    // writeToCache is true and all hashes have already been retrieved - save to the cache
    case Event(
    response: JobSucceededResponse,
    data@ResponsePendingData(_, _, Some(Success(hashes)), _, _, _, _)
    ) if effectiveCallCachingMode.writeToCache =>
      eventList ++= response.executionEvents
      // Publish the image used now that we have it as we might lose the information if Cromwell is restarted
      // in between writing to the cache and writing to the job store
      response.dockerImageUsed foreach publishDockerImageUsed
      saveCacheResults(hashes, data.withSuccessResponse(response))
    // Hashes are still missing and we want them (writeToCache is true) - wait for them
    case Event(response: JobSucceededResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      eventList ++= response.executionEvents
      stay using data.withSuccessResponse(response)
    // Hashes are missing but writeToCache is OFF - complete the job
    case Event(response: JobSucceededResponse, data: ResponsePendingData) =>
      eventList ++= response.executionEvents
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
  }

  // Handles BackendJobFailedResponse messages
  val jobFailedHandler: StateFunction = {
    // writeToCache is true and all hashes already retrieved - save to job store
    case Event(
    response: BackendJobFailedResponse,
    data@ResponsePendingData(_, _, Some(Success(_)), _, _, _, _)
    ) if effectiveCallCachingMode.writeToCache =>
      saveJobCompletionToJobStore(data.withFailedResponse(response))
    // Hashes are still missing and we want them (writeToCache is true) - wait for them
    case Event(response: BackendJobFailedResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      stay using data.withFailedResponse(response)
    // Hashes are missing but writeToCache is OFF - complete the job
    case Event(response: BackendJobFailedResponse, data: ResponsePendingData) =>
      saveJobCompletionToJobStore(data.withFailedResponse(response))
  }

  // Handles JobAbortedResponse messages
  val jobAbortedHandler: StateFunction = {
    // If the job is aborted, don't even wait for the hashes we want to abort it ASAP
    case Event(response: JobAbortedResponse, data: ResponsePendingData) =>
      // If we happen to have hashes it doesn't hurt to publish them
      publishHashesToMetadata(data.hashes)
      // Our work is done, aborted jobs don't get saved to the job store
      forwardAndStop(response)
  }

  // Handles hash success responses
  val hashSuccessResponseHandler: StateFunction = {
    // We're getting the hashes and the job has already completed successfully, save to the cache
    case Event(hashes: CallCacheHashes, data: SucceededResponseData) =>
      saveCacheResults(hashes, data)
    // We're getting the hashes and the job has already completed and failed, save to job store
    case Event(hashes: CallCacheHashes, data: FailedResponseData) =>
      saveJobCompletionToJobStore(data.withHashes(Option(Success(hashes))))
    // We're getting the hashes and the job has already completed and aborted.
    // Since we're not waiting for hashes when the job is aborted this should not happen but if it does,
    // publish the hashes and terminate
    case Event(hashes: CallCacheHashes, data: AbortedResponseData) =>
      publishHashesToMetadata(Option(Success(hashes)))
      // Our work is done, aborted jobs don't get saved to the job store
      forwardAndStop(data.response)
    // We're getting the hashes and the job has not completed yet, save them and stay where we are
    case Event(hashes: CallCacheHashes, data) =>
      addHashesAndStay(data, hashes)
    // Not sure why this is here but pre-existing so leaving it, there's nothing we can do with it anyway at this point
    case Event(CacheMiss, _) =>
      stay()
    case Event(_: CacheHit, _) =>
      stay()
  }

  // Handles hash failure responses
  val hashFailureResponseHandler: StateFunction = {
    // We're getting hash errors but the job was already successful, disable call caching and save to job store
    case Event(HashError(t), data: SucceededResponseData) =>
      disableCallCaching(Option(t))
      saveJobCompletionToJobStore(data.copy(hashes = Option(Failure(t))))
    // We're getting hash errors but the job was already failed, disable call caching and save to job store
    case Event(HashError(t), data: FailedResponseData) =>
      disableCallCaching(Option(t))
      saveJobCompletionToJobStore(data.copy(hashes = Option(Failure(t))))
    // We're getting hash errors but the job was already failed, disable call caching and terminate
    case Event(HashError(t), data: AbortedResponseData) =>
      disableCallCaching(Option(t))
      forwardAndStop(data.response)
    // We're getting hash errors and the job is still running, disable call caching and stay here to wait for the job to finish
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCallCaching(Option(t))
      stay using data.copy(hashes = Option(Failure(t)))
  }

  when(RunningJob)(jobSuccessHandler.orElse(jobFailedHandler).orElse(jobAbortedHandler).orElse(hashSuccessResponseHandler).orElse(hashFailureResponseHandler))

  // When UpdatingCallCache, the FSM always has SucceededResponseData.
  when(UpdatingCallCache) {
    case Event(CallCacheWriteSuccess, data: SucceededResponseData) =>
      saveJobCompletionToJobStore(data)
    case Event(CallCacheWriteFailure(reason), data: SucceededResponseData) =>
      log.error(reason, "{}: Failure writing to call cache: {}", jobTag, reason.getMessage)
      saveJobCompletionToJobStore(data)
  }

  // When UpdatingJobStore, the FSM always has ResponseData.
  when(UpdatingJobStore) {
    case Event(JobStoreWriteSuccess(_), data: ResponseData) =>
      forwardAndStop(data.response)
    case Event(JobStoreWriteFailure(t), _: ResponseData) =>
      respondAndStop(JobFailedNonRetryableResponse(jobDescriptorKey, new Exception(s"JobStore write failure: ${t.getMessage}", t), None))
  }

  onTransition {
    case fromState -> toState =>
      log.debug("Transitioning from {}({}) to {}({})", fromState, stateData, toState, nextStateData)

      EngineJobExecutionActorState.transitionEventString(fromState, toState) foreach {
        eventList :+= ExecutionEvent(_)
      }

  }

  whenUnhandled {
    case Event(EngineStatsActor.JobCountQuery, _) =>
      sender ! EngineStatsActor.JobCount(1)
      stay()
    case Event(e: ActorInitializationException, _) =>
      respondAndStop(JobFailedNonRetryableResponse(jobDescriptorKey, e, None))
    case Event(EngineLifecycleActorAbortCommand, data: ResponsePendingData) =>
      data.backendJobActor match {
        // If there's a backend actor running, send an abort command and wait to hear back
        case Some(backendActor) =>
          backendActor ! AbortJobCommand
          stay()
        // Otherwise we're good (it would be nice to persist the fact that this job is aborted in the job store, but the
        // job store doesn't support that yet)
        case None =>
          forwardAndStop(JobAbortedResponse(jobDescriptorKey))
      }
    case Event(EngineLifecycleActorAbortCommand, _) =>
      forwardAndStop(JobAbortedResponse(jobDescriptorKey))
    case Event(_: CopyingOutputsFailedResponse, _) =>
      // Normally these are caught in when(BackendIsCopyingCachedOutputs) { ... }
      // But sometimes, we get failed responses long after we've moved on from a particular cache hit (usually
      // due to timeouts). That's ok, we just ignore this message in any other situation:
      stay()
    case Event(msg, _) =>
      log.error("Bad message from {} to EngineJobExecutionActor in state {}(with data {}): {}", sender, stateName, stateData, msg)
      stay
  }

  private def publishHashesToMetadata(maybeHashes: Option[Try[CallCacheHashes]]) = publishHashResultsToMetadata(maybeHashes.map(_.map(_.hashes)))
  private def publishDockerImageUsed(image: String) = writeToMetadata(Map("dockerImageUsed" -> image))
  private def publishHashResultsToMetadata(maybeHashes: Option[Try[Set[HashResult]]]) = maybeHashes match {
    case Some(Success(hashes)) =>
      val hashMap = hashes.collect({
        case HashResult(HashKey(useInCallCaching, keyComponents), HashValue(value)) if useInCallCaching =>
          (callCachingHashes + MetadataKey.KeySeparator + keyComponents.mkString(MetadataKey.KeySeparator.toString)) -> value
      }).toMap
      writeToMetadata(hashMap)
    case _ =>
  }

  private def handleReadFromCacheOn(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, updatedData: ResponsePendingData) = {
    jobDescriptor.maybeCallCachingEligible match {
      // If the job is eligible, initialize job hashing and go to CheckingCallCache state
      case eligible: CallCachingEligible =>
        initializeJobHashing(jobDescriptor, activity, eligible) match {
          case Success(ejha) => goto(CheckingCallCache) using updatedData.withEJHA(ejha)
          case Failure(failure) => respondAndStop(JobFailedNonRetryableResponse(jobDescriptorKey, failure, None))
        }
      case _ =>
        // If the job is ineligible, turn call caching off
        writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Cache Miss"))
        disableCallCaching()
        runJob(updatedData)
    }
  }

  private def handleReadFromCacheOff(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, updatedData: ResponsePendingData) = {
    jobDescriptor.maybeCallCachingEligible match {
      // If the job is eligible, initialize job hashing so it can be written to the cache
      case eligible: CallCachingEligible => initializeJobHashing(jobDescriptor, activity, eligible) match {
        case Failure(failure) => log.error(failure, "Failed to initialize job hashing. The job will not be written to the cache")
        case _ =>
      }
      // Don't even initialize hashing to write to the cache if the job is ineligible
      case _ => disableCallCaching()
    }
    // If read from cache is off, always run the job
    runJob(updatedData)
  }

  private def requestExecutionToken(): Unit = {
    jobTokenDispenserActor ! JobExecutionTokenRequest(workflowDescriptor.backendDescriptor.hogGroup, backendLifecycleActorFactory.jobExecutionTokenType)
  }

  // Return the execution token (if we have one)
  private def returnExecutionToken(): Unit = if (stateName != Pending && stateName != RequestingExecutionToken) {
    jobTokenDispenserActor ! JobExecutionTokenReturn
  }

  private def forwardAndStop(response: BackendJobExecutionResponse): State = {
    replyTo forward response
    stop(response)
  }

  private def respondAndStop(response: BackendJobExecutionResponse): State = {
    replyTo ! response
    stop(response)
  }

  private def stop(response: BackendJobExecutionResponse): State = {
    returnExecutionToken()
    instrumentJobComplete(response)
    pushExecutionEventsToMetadataService(jobDescriptorKey, eventList)
    recordExecutionStepTiming(stateName.toString, currentStateDuration)
    context stop self
    stay()
  }

  // Note: StatsD will automatically add a counter value so ne need to separately increment a counter.
  private def instrumentJobComplete(response: BackendJobExecutionResponse) = {
    setJobTimePerState(response, (System.currentTimeMillis() - jobStartTime).millis)
  }

  private def disableCallCaching(reason: Option[Throwable] = None) = {
    reason foreach { e => log.error("{}: Hash error ({}), disabling call caching for this job.", jobTag, e.getMessage) }
    effectiveCallCachingMode = CallCachingOff
    writeCallCachingModeToMetadata()
    writeToMetadata(Map(callCachingHitResultMetadataKey -> false))
  }

  private def disableCacheWrite(reason: Throwable) = {
    log.error(reason, "{}: Disabling cache writing for this job.", jobTag)
    if (effectiveCallCachingMode.writeToCache) {
      effectiveCallCachingMode = effectiveCallCachingMode.withoutWrite
      writeCallCachingModeToMetadata()
    }
  }

  def writeCallCachingModeToMetadata(): Unit = {
    writeToMetadata(Map(effectiveCallCachingKey -> effectiveCallCachingMode.toString))
    writeToMetadata(Map(callCachingAllowReuseMetadataKey -> effectiveCallCachingMode.writeToCache))
  }

  def createJobPreparationActor(jobPrepProps: Props, name: String): ActorRef = context.actorOf(jobPrepProps, name)
  def prepareJob(valueStore: ValueStore) = {
    writeCallCachingModeToMetadata()
    val jobPreparationActorName = s"BackendPreparationActor_for_$jobTag"
    val jobPrepProps = JobPreparationActor.props(workflowDescriptor, jobDescriptorKey, backendLifecycleActorFactory, workflowDockerLookupActor = workflowDockerLookupActor,
      initializationData, serviceRegistryActor = serviceRegistryActor, ioActor = ioActor, backendSingletonActor = backendSingletonActor)
    val jobPreparationActor = createJobPreparationActor(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! CallPreparation.Start(valueStore)
    goto(PreparingJob)
  }

  def requestValueStore() = {
    replyTo ! RequestValueStore
    goto(WaitingForValueStore)
  }

  def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, callCachingEligible: CallCachingEligible): Try[ActorRef] = {
    val maybeFileHashingActorProps = backendLifecycleActorFactory.fileHashingActorProps map {
      _.apply(jobDescriptor, initializationData, serviceRegistryActor, ioActor, fileHashCachingActor)
    }

    maybeFileHashingActorProps match {
      case Some(fileHashingActorProps) =>
        val props = EngineJobHashingActor.props(
          self,
          serviceRegistryActor,
          jobDescriptor,
          initializationData,
          fileHashingActorProps,
          CallCacheReadingJobActor.props(callCacheReadActor, callCachePathPrefixes),
          backendLifecycleActorFactory.runtimeAttributeDefinitions(initializationData),
          backendLifecycleActorFactory.nameForCallCachingPurposes,
          activity,
          callCachingEligible,
          callCachePathPrefixes
        )
        val ejha = context.actorOf(props, s"ejha_for_$jobDescriptor")

        Success(ejha)
      case None => Failure(new IllegalStateException("Tried to initialize job hashing without a file hashing actor !"))
    }
  }

  def makeFetchCachedResultsActor(callCachingEntryId: CallCachingEntryId): Unit = {
    context.actorOf(FetchCachedResultsActor.props(callCachingEntryId, self,
      new CallCache(EngineServicesStore.engineDatabaseInterface)))
    ()
  }

  private def fetchCachedResults( callCachingEntryId: CallCachingEntryId, data: ResponsePendingData) = {
    makeFetchCachedResultsActor(callCachingEntryId)
    goto(FetchingCachedOutputsFromDatabase) using data
  }

  private def makeBackendCopyCacheHit(womValueSimpletons: Seq[WomValueSimpleton],
                                      jobDetritusFiles: Map[String,String],
                                      returnCode: Option[Int],
                                      data: ResponsePendingData,
                                      cacheResultId: CallCachingEntryId,
                                      cacheCopyAttempt: Int) = {
    backendLifecycleActorFactory.cacheHitCopyingActorProps match {
      case Some(propsMaker) =>
        val backendCacheHitCopyingActorProps = propsMaker(data.jobDescriptor, initializationData, serviceRegistryActor, ioActor, cacheCopyAttempt, blacklistCache)
        val cacheHitCopyActor = context.actorOf(backendCacheHitCopyingActorProps, buildCacheHitCopyingActorName(data.jobDescriptor, cacheResultId))
        cacheHitCopyActor ! CopyOutputsCommand(womValueSimpletons, jobDetritusFiles, returnCode)
        replyTo ! JobRunning(data.jobDescriptor.key, data.jobDescriptor.evaluatedTaskInputs)
        goto(BackendIsCopyingCachedOutputs)
      case None =>
        // This should be impossible with the FSM, but luckily, we CAN recover if some foolish future programmer makes this happen:
        val errorMessage = "Call caching copying should never have even been attempted with no copy actor props! (Programmer error!)"
        log.error(errorMessage)
        self ! JobFailedNonRetryableResponse(data.jobDescriptor.key, new RuntimeException(errorMessage), None)
        goto(BackendIsCopyingCachedOutputs)
    }
  }

  private [job] def createBackendJobExecutionActor(data: ResponsePendingData) = {
    context.actorOf(data.bjeaProps, BackendJobExecutionActor.buildJobExecutionActorName(workflowIdForLogging, data.jobDescriptor.key))
  }

  private def runJob(data: ResponsePendingData) = {
    val backendJobExecutionActor = createBackendJobExecutionActor(data)
    backendJobExecutionActor ! command
    replyTo ! JobRunning(data.jobDescriptor.key, data.jobDescriptor.evaluatedTaskInputs)
    goto(RunningJob) using data.withBackendActor(backendJobExecutionActor)
  }

  private def handleCacheInvalidatedResponse(response: CallCacheInvalidatedResponse, data: ResponsePendingData) = {
    def updateMetadataForInvalidatedEntry(entry: CallCachingEntry) = {
      import cromwell.core.ExecutionIndex._
      import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter

      val workflowId = WorkflowId.fromString(entry.workflowExecutionUuid)
      // If the entry doesn't have an attempt, it means that this cache entry was added before this change
      // and we don't know which attempt yielded this cache entry
      // In that case make a best effort and update the first attempt
      val key = Option((entry.callFullyQualifiedName, entry.jobIndex.toIndex, entry.jobAttempt.getOrElse(1)))
      serviceRegistryActor.putMetadataWithRawKey(workflowId, key, Map(callCachingAllowReuseMetadataKey -> false))
    }

    response match {
      case CallCacheInvalidatedFailure(failure) => log.error(failure, "Failed to invalidate cache entry for job: {}", jobDescriptorKey)
      case CallCacheInvalidatedSuccess(Some(entry)) => updateMetadataForInvalidatedEntry(entry)
      case _ =>
    }

    data.ejha match {
      case Some(ejha) =>
        workflowLogger.debug("Trying to use another cache hit for job: {}", jobDescriptorKey)
        ejha ! NextHit
        goto(CheckingCallCache)
      case _ =>
        workflowLogger.info(
          "Could not find a suitable cache hit. " +
            "Call cache hit process had {} total hit failures before completing unsuccessfully. " +
            "Falling back to running job: {}", data.cacheHitFailureCount, jobDescriptorKey)
        runJob(data)
    }
  }

  private def buildCacheHitCopyingActorName(jobDescriptor: BackendJobDescriptor, cacheResultId: CallCachingEntryId) = {
    s"$workflowIdForLogging-BackendCacheHitCopyingActor-$jobTag-${cacheResultId.id}"
  }

  private def logCacheHitSuccess(data: ResponsePendingData): Unit = {
    workflowLogger.info(
      "Call cache hit process had {} total hit failures before completing successfully",
      data.cacheHitFailureCount,
    )
  }

  private def logCacheHitFailure(data: ResponsePendingData, reason: Throwable): Unit = {
    val problemSummary =
      s"Failed copying cache results for job $jobDescriptorKey (${reason.getClass.getSimpleName}: ${reason.getMessage})"
    if (invalidationRequired) {
      // Whenever invalidating a cache result, always log why the invalidation occurred
      workflowLogger.warn(s"$problemSummary, invalidating cache entry.")
    } else if (data.cacheHitFailureCount < 3) {
      workflowLogger.info(problemSummary)
    }
  }

  private def invalidateCacheHitAndTransition(ejeaCacheHit: EJEACacheHit, data: ResponsePendingData, reason: Throwable) = {
    logCacheHitFailure(data, reason)
    val updatedData = data.copy(cacheHitFailureCount = data.cacheHitFailureCount + 1)

    if (invalidationRequired) {
      invalidateCacheHit(ejeaCacheHit.hit.cacheResultId)
      goto(InvalidatingCacheEntry) using updatedData
    } else {
      handleCacheInvalidatedResponse(CallCacheInvalidationUnnecessary, updatedData)
    }
  }

  protected def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = {
    val callCache = new CallCache(EngineServicesStore.engineDatabaseInterface)
    context.actorOf(CallCacheInvalidateActor.props(callCache, cacheId), s"CallCacheInvalidateActor${cacheId.id}-$tag")
    ()
  }

  private def checkCacheEntryExistence() = {
    callCacheReadActor ! CallCacheEntryForCall(workflowIdForLogging, jobDescriptorKey)
    goto(CheckingCacheEntryExistence)
  }

  private def saveCacheResults(hashes: CallCacheHashes, data: SucceededResponseData) = {
    callCacheWriteActor ! SaveCallCacheHashes(CallCacheHashBundle(workflowIdForLogging, hashes, data.response))
    val updatedData = data.copy(hashes = Option(Success(hashes)))
    goto(UpdatingCallCache) using updatedData
  }

  private def saveJobCompletionToJobStore(updatedData: ShouldBeSavedToJobStoreResponseData) = {
    updatedData match {
      case SucceededResponseData(JobSucceededResponse(jobKey, returnCode, jobOutputs, _, _, _, _), hashes) =>
        publishHashesToMetadata(hashes)
        saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case FailedResponseData(JobFailedNonRetryableResponse(jobKey, throwable, returnCode), hashes) =>
        publishHashesToMetadata(hashes)
        writeToMetadata(Map(callCachingAllowReuseMetadataKey -> false))
        saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case FailedResponseData(JobFailedRetryableResponse(jobKey, throwable, returnCode), hashes) =>
        publishHashesToMetadata(hashes)
        writeToMetadata(Map(callCachingAllowReuseMetadataKey -> false))
        saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }

    updatedData.dockerImageUsed foreach publishDockerImageUsed
    goto(UpdatingJobStore) using updatedData
  }

  private def saveSuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], outputs: CallOutputs) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowIdForLogging)
    val jobStoreResult = JobResultSuccess(returnCode, outputs)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }

  private def saveUnsuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], reason: Throwable, retryable: Boolean) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowIdForLogging)
    val jobStoreResult = JobResultFailure(returnCode, reason, retryable)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }

  private def writeToMetadata(keyValues: Map[String, Any]) = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(workflowIdForLogging, Option(jobDescriptorKey), keyValues)
  }

  private def addHashesAndStay(data: EJEAData, hashes: CallCacheHashes): State = {
    val updatedData = data match {
      case responsePending: ResponsePendingData => responsePending.copy(hashes = Option(Success(hashes)))
      case responseData: ResponseData => responseData.withHashes(Option(Success(hashes)))
      case _ => data
    }
    stay using updatedData
  }
}

object EngineJobExecutionActor {
  /** States */
  sealed trait EngineJobExecutionActorState
  case object Pending extends EngineJobExecutionActorState
  case object RequestingExecutionToken extends EngineJobExecutionActorState
  case object CheckingJobStore extends EngineJobExecutionActorState
  case object CheckingCallCache extends EngineJobExecutionActorState
  case object FetchingCachedOutputsFromDatabase extends EngineJobExecutionActorState
  case object BackendIsCopyingCachedOutputs extends EngineJobExecutionActorState
  case object WaitingForValueStore extends EngineJobExecutionActorState
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState
  case object UpdatingCallCache extends EngineJobExecutionActorState
  case object CheckingCacheEntryExistence extends EngineJobExecutionActorState
  case object UpdatingJobStore extends EngineJobExecutionActorState
  case object InvalidatingCacheEntry extends EngineJobExecutionActorState

  object EngineJobExecutionActorState {
    def transitionEventString(fromState: EngineJobExecutionActorState, toState: EngineJobExecutionActorState): Option[String] = {

      def callCacheStateGroup: Set[EngineJobExecutionActorState] = Set(
        CheckingCallCache,
        FetchingCachedOutputsFromDatabase,
        BackendIsCopyingCachedOutputs,
        CheckingCacheEntryExistence,
        InvalidatingCacheEntry
      )

      if (fromState == toState) None
      else if (callCacheStateGroup.contains(fromState) && callCacheStateGroup.contains(toState)) None
      else if (callCacheStateGroup.contains(toState)) Option("CallCacheReading")
      else Option(toState.toString)
    }
  }



  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  val CacheMetadataKeyPrefix = CallMetadataKeys.CallCaching + MetadataKey.KeySeparator

  def props(replyTo: ActorRef,
            jobDescriptorKey: BackendJobDescriptorKey,
            workflowDescriptor: EngineWorkflowDescriptor,
            backendLifecycleActorFactory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            restarting: Boolean,
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef,
            jobStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            callCacheWriteActor: ActorRef,
            workflowDockerLookupActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonActor: Option[ActorRef],
            callCachingMode: CallCachingMode,
            command: BackendJobExecutionActorCommand,
            fileHashCacheActor: Option[ActorRef],
            blacklistCache: Option[BlacklistCache]) = {
    Props(new EngineJobExecutionActor(
      replyTo = replyTo,
      jobDescriptorKey = jobDescriptorKey,
      workflowDescriptor = workflowDescriptor,
      backendLifecycleActorFactory = backendLifecycleActorFactory,
      initializationData = initializationData,
      restarting = restarting,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      jobStoreActor = jobStoreActor,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      workflowDockerLookupActor = workflowDockerLookupActor,
      jobTokenDispenserActor = jobTokenDispenserActor,
      backendSingletonActor = backendSingletonActor,
      callCachingMode = callCachingMode,
      command = command,
      fileHashCachingActor = fileHashCacheActor,
      blacklistCache = blacklistCache)).withDispatcher(EngineDispatcher)
  }

  case class EJEACacheHit(hit: CacheHit, hitNumber: Int, details: Option[String])

  private[execution] sealed trait EJEAData {
    override def toString = getClass.getSimpleName
  }

  private[execution] case object NoData extends EJEAData

  private[execution] case class ResponsePendingData(jobDescriptor: BackendJobDescriptor,
                                                    bjeaProps: Props,
                                                    hashes: Option[Try[CallCacheHashes]] = None,
                                                    ejha: Option[ActorRef] = None,
                                                    ejeaCacheHit: Option[EJEACacheHit] = None,
                                                    backendJobActor: Option[ActorRef] = None,
                                                    cacheHitFailureCount: Int = 0
                                                   ) extends EJEAData {

    def withEJHA(ejha: ActorRef): EJEAData = this.copy(ejha = Option(ejha))

    def withBackendActor(actorRef: ActorRef) = this.copy(backendJobActor = Option(actorRef))

    def withSuccessResponse(success: JobSucceededResponse): SucceededResponseData = SucceededResponseData(success, hashes)
    def withFailedResponse(failed: BackendJobFailedResponse): FailedResponseData = FailedResponseData(failed, hashes)
    def withAbortedResponse(aborted: JobAbortedResponse): AbortedResponseData = AbortedResponseData(aborted, hashes)

    def withCacheHit(cacheHit: CacheHit) = {
      val newEjeaCacheHit = ejeaCacheHit map { currentCacheHit =>
        currentCacheHit.copy(hit = cacheHit, hitNumber = currentCacheHit.hitNumber + 1)
      } getOrElse EJEACacheHit(cacheHit, 0, None)

      this.copy(ejeaCacheHit = Option(newEjeaCacheHit))
    }

    def withCacheDetails(details: String) = this.copy(ejeaCacheHit = ejeaCacheHit.map(_.copy(details = Option(details))))
  }

  private[execution] trait ResponseData extends EJEAData {
    def response: BackendJobExecutionResponse
    def hashes: Option[Try[CallCacheHashes]]
    def dockerImageUsed: Option[String]
    def withHashes(hashes: Option[Try[CallCacheHashes]]): ResponseData
  }

  // Only Successes and Failures are saved to the job store, not Aborts. Why ? Because. This could be an improvement AFAICT.
  private[execution] trait ShouldBeSavedToJobStoreResponseData extends ResponseData

  private[execution] case class SucceededResponseData(successResponse: JobSucceededResponse,
                                                       hashes: Option[Try[CallCacheHashes]] = None) extends ShouldBeSavedToJobStoreResponseData {
    override def response = successResponse
    override def dockerImageUsed = successResponse.dockerImageUsed
    override def withHashes(hashes: Option[Try[CallCacheHashes]]) = copy(hashes = hashes)
  }

  private[execution] case class FailedResponseData(failedResponse: BackendJobFailedResponse,
                                                      hashes: Option[Try[CallCacheHashes]] = None) extends ShouldBeSavedToJobStoreResponseData {
    override def response = failedResponse
    // Seems like we should be able to get the docker image used even if the job failed
    override def dockerImageUsed = None
    override def withHashes(hashes: Option[Try[CallCacheHashes]]) = copy(hashes = hashes)
  }

  private[execution] case class AbortedResponseData(abortedResponse: JobAbortedResponse,
                                                   hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData {
    override def response = abortedResponse
    override def dockerImageUsed = None
    override def withHashes(hashes: Option[Try[CallCacheHashes]]) = copy(hashes = hashes)
  }
}
