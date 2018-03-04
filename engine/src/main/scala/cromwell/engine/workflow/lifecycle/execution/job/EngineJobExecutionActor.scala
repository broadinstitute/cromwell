package cromwell.engine.workflow.lifecycle.execution.job

import java.util.concurrent.TimeoutException

import akka.actor.SupervisorStrategy.{Escalate, Stop}
import akka.actor.{ActorInitializationException, ActorRef, LoggingFSM, OneForOneStrategy, Props}
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.database.sql.tables.CallCachingEntry
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.instrumentation.JobInstrumentation
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestValueStore
import cromwell.engine.workflow.lifecycle.execution._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCacheHashBundle
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CallCacheHashes, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed}
import cromwell.engine.workflow.lifecycle.execution.job.preparation.{CallPreparation, JobPreparationActor}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore._
import cromwell.services.EngineServicesStore
import cromwell.services.metadata.CallMetadataKeys.CallCachingKeys
import cromwell.services.metadata.{CallMetadataKeys, MetadataJobKey, MetadataKey}
import cromwell.webservice.EngineStatsActor

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

class EngineJobExecutionActor(replyTo: ActorRef,
                              jobDescriptorKey: BackendJobDescriptorKey,
                              workflowDescriptor: EngineWorkflowDescriptor,
                              factory: BackendLifecycleActorFactory,
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
                              backendName: String,
                              callCachingMode: CallCachingMode,
                              command: BackendJobExecutionActorCommand) extends LoggingFSM[EngineJobExecutionActorState, EJEAData]
  with WorkflowLogging with CallMetadataHelper with JobInstrumentation {

  override val workflowIdForLogging = workflowDescriptor.id
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

  // There's no need to check for a cache hit again if we got preempted, or if there's no result copying actor defined
  // NB: this can also change (e.g. if we have a HashError we just force this to CallCachingOff)
  private[execution] var effectiveCallCachingMode = {
    if (factory.fileHashingActorProps.isEmpty) CallCachingOff
    else if (factory.cacheHitCopyingActorProps.isEmpty || jobDescriptorKey.attempt > 1) {
      callCachingMode.withoutRead
    } else callCachingMode
  }

  // For tests:
  private[execution] def checkEffectiveCallCachingMode = effectiveCallCachingMode

  private val effectiveCallCachingKey = CallCachingKeys.EffectiveModeKey
  private val callCachingReadResultMetadataKey = CallCachingKeys.ReadResultMetadataKey
  private val callCachingHitResultMetadataKey = CallCachingKeys.HitResultMetadataKey
  private val callCachingAllowReuseMetadataKey = CallCachingKeys.AllowReuseMetadataKey
  private val callCachingHitFailures = CallCachingKeys.HitFailuresKey
  private val callCachingHashes = CallCachingKeys.HashesKey

  implicit val ec: ExecutionContext = context.dispatcher

  override def preStart() = {
    log.debug(s"$tag: $effectiveCallCachingKey: $effectiveCallCachingMode")
  }

  startWith(Pending, NoData)
  private var eventList: Seq[ExecutionEvent] = Seq(ExecutionEvent(stateName.toString))

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
  // writing to the JobStore. In that case, we don't want to cache to ourselves (turn cache read off), nor do we want to 
  // try and write the cache info again, which would fail (turn cache write off).
  // This means call caching should be disabled.
  // Note that we check that there is not a Cache entry for *this* current job. It's still technically possible
  // to call cache to another job that finished while this one was running (before the restart).
  when(CheckingCacheEntryExistence) {
    // There was already a cache entry for this job
    case Event(HasCallCacheEntry(_), NoData) =>
      // Disable call caching
      effectiveCallCachingMode = CallCachingOff
      requestValueStore()
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
    case Event(CachedOutputLookupSucceeded(womValueSimpletons, jobDetritus, returnCode, cacheResultId, cacheHitDetails), data: ResponsePendingData) =>
      writeToMetadata(Map(
        callCachingHitResultMetadataKey -> true,
        callCachingReadResultMetadataKey -> s"Cache Hit: $cacheHitDetails"))
      log.debug("Cache hit for {}! Fetching cached result {}", jobTag, cacheResultId)
      makeBackendCopyCacheHit(womValueSimpletons, jobDetritus, returnCode, data, cacheResultId) using data.withCacheDetails(cacheHitDetails)
    case Event(CachedOutputLookupFailed(_, error), data: ResponsePendingData) =>
      log.warning("Can't make a copy of the cached job outputs for {} due to {}. Running job.", jobTag, error)
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
    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)), _, _, _)) =>
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: JobSucceededResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      // Wait for the CallCacheHashes
      stay using data.withSuccessResponse(response)
    case Event(response: JobSucceededResponse, data: ResponsePendingData) => // bad hashes or cache write off
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
    case Event(response: BackendJobExecutionResponse, data @ ResponsePendingData(_, _, _, _, Some(cacheHit), _)) =>
      response match {
        case f: BackendJobFailedResponse => invalidateCacheHitAndTransition(cacheHit, data, f.throwable)
        case _ => runJob(data)
      }

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

  when(RunningJob) {
    case Event(hashes: CallCacheHashes, data: SucceededResponseData) =>
      saveCacheResults(hashes, data)
    case Event(hashes: CallCacheHashes, data) =>
      addHashesAndStay(data, hashes)
    case Event(CacheMiss, _) =>
      stay()
    case Event(_: CacheHit, _) =>
      stay()

    case Event(HashError(t), data: SucceededResponseData) =>
      disableCallCaching(Option(t))
      saveJobCompletionToJobStore(data.copy(hashes = Option(Failure(t))))
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCallCaching(Option(t))
      stay using data.copy(hashes = Option(Failure(t)))

    // Success
    // All hashes already retrieved - save to the cache
    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)), _, _, _)) if effectiveCallCachingMode.writeToCache =>
      eventList ++= response.executionEvents
      saveCacheResults(hashes, data.withSuccessResponse(response))
    // Some hashes are still missing - waiting for them
    case Event(response: JobSucceededResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      eventList ++= response.executionEvents
      log.debug(s"Got job result for {}, awaiting hashes", jobTag)
      stay using data.withSuccessResponse(response)
    // writeToCache is OFF - complete the job
    case Event(response: JobSucceededResponse, data: ResponsePendingData) =>
      eventList ++= response.executionEvents
      saveJobCompletionToJobStore(data.withSuccessResponse(response))

    // Non-success:
    // Even if the job failed, wait for the hashes so we can publish them to metadata
    case Event(response: BackendJobExecutionResponse, data: ResponsePendingData) if effectiveCallCachingMode.writeToCache && data.hashes.isEmpty =>
      stay() using data.withResponse(response)
    case Event(response: BackendJobExecutionResponse, data: ResponsePendingData) =>
      saveJobCompletionToJobStore(data.withResponse(response))
  }

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
      eventList :+= ExecutionEvent(toState.toString)
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
    // We can get extra I/O timeout messages from previous cache copy attempt. This is due to the fact that if one file fails to copy,
    // we immediately cache miss and try the next potential hit, but don't cancel previous copy requests.
    case Event(JobFailedNonRetryableResponse(_, _: TimeoutException, _), _) =>
      stay()
    case Event(msg, _) =>
      log.error("Bad message from {} to EngineJobExecutionActor in state {}(with data {}): {}", sender, stateName, stateData, msg)
      stay
  }

  private def publishHashesToMetadata(maybeHashes: Option[Try[CallCacheHashes]]) = maybeHashes match {
    case Some(Success(hashes)) =>
      val hashMap = hashes.hashes.collect({
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
    jobTokenDispenserActor ! JobExecutionTokenRequest(factory.jobExecutionTokenType)
  }

  private def returnExecutionToken(): Unit = {
    jobTokenDispenserActor ! JobExecutionTokenReturn
  }

  private def forwardAndStop(response: BackendJobExecutionResponse): State = {
    replyTo forward response
    returnExecutionToken()
    instrumentJobComplete(response)
    pushExecutionEventsToMetadataService(jobDescriptorKey, eventList)
    context stop self
    stay()
  }

  private def respondAndStop(response: BackendJobExecutionResponse): State = {
    replyTo ! response
    returnExecutionToken()
    instrumentJobComplete(response)
    pushExecutionEventsToMetadataService(jobDescriptorKey, eventList)
    context stop self
    stay()
  }
  
  // Note: StatsD will automatically add a counter value so ne need to separately increment a counter.
  private def instrumentJobComplete(response: BackendJobExecutionResponse) = {
    setJobTimePerState(response, (System.currentTimeMillis() - jobStartTime).millis)
  }

  private def disableCallCaching(reason: Option[Throwable] = None) = {
    reason foreach { log.error(_, "{}: Hash error, disabling call caching for this job.", jobTag) }
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
    val jobPrepProps = JobPreparationActor.props(workflowDescriptor, jobDescriptorKey, factory, workflowDockerLookupActor = workflowDockerLookupActor,
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
    val maybeFileHashingActorProps = factory.fileHashingActorProps map {
      _.apply(jobDescriptor, initializationData, serviceRegistryActor, ioActor)
    }

    maybeFileHashingActorProps match {
      case Some(fileHashingActorProps) =>
        val props = EngineJobHashingActor.props(
          self,
          serviceRegistryActor,
          jobDescriptor,
          initializationData,
          fileHashingActorProps,
          CallCacheReadingJobActor.props(callCacheReadActor),
          factory.runtimeAttributeDefinitions(initializationData),
          backendName,
          activity,
          callCachingEligible
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

  private def makeBackendCopyCacheHit(womValueSimpletons: Seq[WomValueSimpleton], jobDetritusFiles: Map[String,String], returnCode: Option[Int], data: ResponsePendingData, cacheResultId: CallCachingEntryId) = {
    factory.cacheHitCopyingActorProps match {
      case Some(propsMaker) =>
        val backendCacheHitCopyingActorProps = propsMaker(data.jobDescriptor, initializationData, serviceRegistryActor, ioActor)
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
        workflowLogger.info("Could not find a suitable cache hit, falling back to running job: {}", jobDescriptorKey)
        runJob(data)
    }
  }

  private def buildCacheHitCopyingActorName(jobDescriptor: BackendJobDescriptor, cacheResultId: CallCachingEntryId) = {
    s"$workflowIdForLogging-BackendCacheHitCopyingActor-$jobTag-${cacheResultId.id}"
  }

  private def publishHitFailure(cache: EJEACacheHit, failure: Throwable) = {
    import MetadataKey._
    import cromwell.services.metadata.MetadataService._

    cache.details foreach { details =>
      val metadataKey = MetadataKey(
        workflowIdForLogging,
        Option(MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)),
        s"$callCachingHitFailures[${cache.hitNumber}]:${details.escapeMeta}"
      )

      serviceRegistryActor ! PutMetadataAction(throwableToMetadataEvents(metadataKey, failure))
    }
  }

  private def invalidateCacheHitAndTransition(ejeaCacheHit: EJEACacheHit, data: ResponsePendingData, reason: Throwable) = {
    publishHitFailure(ejeaCacheHit, reason)

    val invalidationRequired = effectiveCallCachingMode match {
      case CallCachingOff => throw new RuntimeException("Should not be calling invalidateCacheHit if call caching is off!") // Very unexpected. Fail out of this bad-state EJEA.
      case activity: CallCachingActivity => activity.options.invalidateBadCacheResults
    }
    if (invalidationRequired) {
      log.error(reason, "Failed copying cache results for job {}, invalidating cache entry.", jobDescriptorKey)
      invalidateCacheHit(ejeaCacheHit.hit.cacheResultId)
      goto(InvalidatingCacheEntry)
    } else {
      handleCacheInvalidatedResponse(CallCacheInvalidationUnnecessary, data)
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

  private def saveJobCompletionToJobStore(updatedData: ResponseData) = {
    updatedData.response match {
      case JobSucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: CallOutputs, _, _, _) =>
        publishHashesToMetadata(updatedData.hashes)
        saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case JobAbortedResponse(_: BackendJobDescriptorKey) =>
        log.debug("{}: Won't save aborted job response to JobStore", jobTag)
        forwardAndStop(updatedData.response)
      case JobFailedNonRetryableResponse(jobKey, throwable: Throwable, returnCode: Option[Int]) =>
        publishHashesToMetadata(updatedData.hashes)
        writeToMetadata(Map(callCachingAllowReuseMetadataKey -> false))
        saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case JobFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) =>
        writeToMetadata(Map(callCachingAllowReuseMetadataKey -> false))
        saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }
    updatedData.dockerImageUsed foreach { image => writeToMetadata(Map("dockerImageUsed" -> image)) }
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

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  val CacheMetadataKeyPrefix = CallMetadataKeys.CallCaching + MetadataKey.KeySeparator

  def props(replyTo: ActorRef,
            jobDescriptorKey: BackendJobDescriptorKey,
            workflowDescriptor: EngineWorkflowDescriptor,
            factory: BackendLifecycleActorFactory,
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
            backendName: String,
            callCachingMode: CallCachingMode,
            command: BackendJobExecutionActorCommand) = {
    Props(new EngineJobExecutionActor(
      replyTo = replyTo,
      jobDescriptorKey = jobDescriptorKey,
      workflowDescriptor = workflowDescriptor,
      factory = factory,
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
      backendName = backendName: String,
      callCachingMode = callCachingMode,
      command = command)).withDispatcher(EngineDispatcher)
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
                                                    backendJobActor: Option[ActorRef] = None
                                                   ) extends EJEAData {

    def withEJHA(ejha: ActorRef): EJEAData = this.copy(ejha = Option(ejha))

    def withBackendActor(actorRef: ActorRef) = this.copy(backendJobActor = Option(actorRef))

    def withSuccessResponse(success: JobSucceededResponse) = SucceededResponseData(success, hashes)

    def withResponse(response: BackendJobExecutionResponse) = response match {
      case success: JobSucceededResponse => SucceededResponseData(success, hashes)
      case failure => NotSucceededResponseData(failure, hashes)
    }

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

  private[execution] case class SucceededResponseData(successResponse: JobSucceededResponse,
                                                      hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData {
    override def response = successResponse
    override def dockerImageUsed = successResponse.dockerImageUsed
    override def withHashes(hashes: Option[Try[CallCacheHashes]]) = copy(hashes = hashes)
  }

  private[execution] case class NotSucceededResponseData(response: BackendJobExecutionResponse,
                                                         hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData {
    override def dockerImageUsed = None
    override def withHashes(hashes: Option[Try[CallCacheHashes]]) = copy(hashes = hashes)
  }
}
