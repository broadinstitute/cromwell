package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql.tables.CallCachingEntry
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CallCacheHashes, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed}
import cromwell.engine.workflow.lifecycle.execution.preparation.{CallPreparation, JobPreparationActor}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDenied, JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore._
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.{CallMetadataKeys, MetadataKey}
import wdl4s.TaskOutput

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success, Try}

class EngineJobExecutionActor(replyTo: ActorRef,
                              jobDescriptorKey: BackendJobDescriptorKey,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              val serviceRegistryActor: ActorRef,
                              ioActor: ActorRef,
                              jobStoreActor: ActorRef,
                              callCacheReadActor: ActorRef,
                              callCacheWriteActor: ActorRef,
                              dockerHashActor: ActorRef,
                              jobTokenDispenserActor: ActorRef,
                              backendSingletonActor: Option[ActorRef],
                              backendName: String,
                              callCachingMode: CallCachingMode) extends LoggingFSM[EngineJobExecutionActorState, EJEAData] with WorkflowLogging with CallMetadataHelper {

  override val workflowIdForLogging = executionData.workflowDescriptor.id
  override val workflowIdForCallMetadata = executionData.workflowDescriptor.id

  val jobTag = s"${workflowIdForLogging.shortString}:${jobDescriptorKey.call.fullyQualifiedName}:${jobDescriptorKey.index.fromIndex}:${jobDescriptorKey.attempt}"
  val tag = s"EJEA_$jobTag"

  // There's no need to check for a cache hit again if we got preempted, or if there's no result copying actor defined
  // NB: this can also change (e.g. if we have a HashError we just force this to CallCachingOff)
  private var effectiveCallCachingMode = {
    if (factory.fileHashingActorProps.isEmpty) CallCachingOff
    else if (factory.cacheHitCopyingActorProps.isEmpty || jobDescriptorKey.attempt > 1) {
      callCachingMode.withoutRead
    } else callCachingMode
  }

  // For tests:
  private[execution] def checkEffectiveCallCachingMode = effectiveCallCachingMode

  private[execution] var executionToken: Option[JobExecutionToken] = None

  private val effectiveCallCachingKey = CacheMetadataKeyPrefix + "effectiveCallCachingMode"
  private val callCachingReadResultMetadataKey = CacheMetadataKeyPrefix + "result"
  private val callCachingHitResultMetadataKey = CacheMetadataKeyPrefix + "hit"
  private val callCachingAllowReuseMetadataKey = CacheMetadataKeyPrefix + "allowResultReuse"

  implicit val ec: ExecutionContext = context.dispatcher

  override def preStart() = {
    log.debug(s"$tag: $effectiveCallCachingKey: $effectiveCallCachingMode")
    writeCallCachingModeToMetadata()
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
    case Event(JobExecutionTokenDispensed(jobExecutionToken), NoData) =>
      executionToken = Option(jobExecutionToken)
      replyTo ! JobStarting(jobDescriptorKey)
      if (restarting) {
        val jobStoreKey = jobDescriptorKey.toJobStoreKey(workflowIdForLogging)
        jobStoreActor ! QueryJobCompletion(jobStoreKey, jobDescriptorKey.call.task.outputs)
        goto(CheckingJobStore)
      } else {
        prepareJob()
      }
    case Event(JobExecutionTokenDenied(positionInQueue), NoData) =>
      log.debug("Token denied so cannot start yet. Currently position {} in the queue", positionInQueue)
      stay()
  }

  // When CheckingJobStore, the FSM always has NoData
  when(CheckingJobStore) {
    case Event(JobNotComplete, NoData) =>
      prepareJob()
    case Event(JobComplete(jobResult), NoData) =>
      val response = jobResult match {
        case JobResultSuccess(returnCode, jobOutputs) => JobSucceededResponse(jobDescriptorKey, returnCode, jobOutputs, None, Seq.empty)
        case JobResultFailure(returnCode, reason, false) => JobFailedNonRetryableResponse(jobDescriptorKey, reason, returnCode)
        case JobResultFailure(returnCode, reason, true) => JobFailedRetryableResponse(jobDescriptorKey, reason, returnCode)
      }
      respondAndStop(response)
    case Event(f: JobStoreReadFailure, NoData) =>
      log.error(f.reason, "{}: Error reading from JobStore", tag)
      // Escalate
      throw new RuntimeException(f.reason)
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
      fetchCachedResults(jobDescriptorKey.call.task.outputs, hit.cacheResultId, data.withCacheHit(hit))
    case Event(HashError(t), data: ResponsePendingData) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Hashing Error: ${t.getMessage}"))
      disableCallCaching(Option(t))
      runJob(data)
  }

  when(FetchingCachedOutputsFromDatabase) {
    case Event(CachedOutputLookupSucceeded(wdlValueSimpletons, jobDetritus, returnCode, cacheResultId, cacheHitDetails), data: ResponsePendingData) =>
      writeToMetadata(Map(
        callCachingHitResultMetadataKey -> true,
        callCachingReadResultMetadataKey -> s"Cache Hit: $cacheHitDetails"))
      log.debug("Cache hit for {}! Fetching cached result {}", jobTag, cacheResultId)
      makeBackendCopyCacheHit(wdlValueSimpletons, jobDetritus, returnCode, data, cacheResultId)
    case Event(CachedOutputLookupFailed(callCachingEntryId, error), data: ResponsePendingData) =>
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
    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)), Some(ejha), _)) =>
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, None, Some(ejha), _)) if effectiveCallCachingMode.writeToCache =>
      // Wait for the CallCacheHashes
      stay using data.withSuccessResponse(response)
    case Event(response: JobSucceededResponse, data: ResponsePendingData) => // bad hashes or cache write off
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
    case Event(response: BackendJobExecutionResponse, data @ ResponsePendingData(_, _, _, _, Some(cacheHit))) =>
      response match {
        case f: BackendJobFailedResponse => invalidateCacheHitAndTransition(cacheHit.cacheResultId, data, f.throwable)
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
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
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

    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)), _, _)) if effectiveCallCachingMode.writeToCache =>
      eventList ++= response.executionEvents
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: JobSucceededResponse, data @ ResponsePendingData(_, _, None, _, _)) if effectiveCallCachingMode.writeToCache =>
      eventList ++= response.executionEvents
      log.debug(s"Got job result for {}, awaiting hashes", jobTag)
      stay using data.withSuccessResponse(response)
    case Event(response: JobSucceededResponse, data: ResponsePendingData) =>
      eventList ++= response.executionEvents
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
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
    case Event(JobStoreWriteFailure(t), data: ResponseData) =>
      respondAndStop(JobFailedNonRetryableResponse(jobDescriptorKey, new Exception(s"JobStore write failure: ${t.getMessage}", t), None))
  }

  onTransition {
    case fromState -> toState =>
      log.debug("Transitioning from {}({}) to {}({})", fromState, stateData, toState, nextStateData)
      eventList :+= ExecutionEvent(toState.toString)
  }

  whenUnhandled {
    case Event(msg, _) =>
      log.error("Bad message from {} to EngineJobExecutionActor in state {}(with data {}): {}", sender, stateName, stateData, msg)
      stay
  }
  
  private def handleReadFromCacheOn(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, updatedData: ResponsePendingData) = {
    jobDescriptor.callCachingEligibility match {
        // If the job is eligible, initialize job hashing and go to CheckingCallCache state
      case CallCachingEligible =>
        initializeJobHashing(jobDescriptor, activity) match {
          case Success(ejha) => goto(CheckingCallCache) using updatedData.withEJHA(ejha)
          case Failure(failure) => respondAndStop(JobFailedNonRetryableResponse(jobDescriptorKey.jobKey, failure, None))
        }
      case ineligible: CallCachingIneligible =>
        // If the job is ineligible, turn call caching off
        writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Cache Miss: ${ineligible.message}"))
        disableCallCaching()
        runJob(updatedData)
    }
  }

  private def handleReadFromCacheOff(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, updatedData: ResponsePendingData) = {
    jobDescriptor.callCachingEligibility match {
        // If the job is eligible, initialize job hashing so it can be written to the cache
      case CallCachingEligible => initializeJobHashing(jobDescriptor, activity) match {
        case Failure(failure) => log.error(failure, "Failed to initialize job hashing. The job will not be written to the cache")
        case _ =>
      }
      // Don't even initialize hashing to write to the cache if the job is ineligible
      case ineligible: CallCachingIneligible => disableCallCaching()
    }
    // If read from cache is off, always run the job
    runJob(updatedData)
  }

  private def requestExecutionToken(): Unit = {
    jobTokenDispenserActor ! JobExecutionTokenRequest(factory.jobExecutionTokenType)
  }

  private def returnExecutionToken(): Unit = {
    executionToken foreach { jobTokenDispenserActor ! JobExecutionTokenReturn(_) }
  }

  private def forwardAndStop(response: Any): State = {
    replyTo forward response
    returnExecutionToken()
    pushExecutionEventsToMetadataService(jobDescriptorKey, eventList)
    context stop self
    stay()
  }

  private def respondAndStop(response: Any): State = {
    replyTo ! response
    returnExecutionToken()
    pushExecutionEventsToMetadataService(jobDescriptorKey, eventList)
    context stop self
    stay()
  }

  private def disableCallCaching(reason: Option[Throwable] = None) = {
    reason foreach { r => log.error(r, "{}: Hash error, disabling call caching for this job.", jobTag) }
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
  def prepareJob() = {
    val jobPreparationActorName = s"BackendPreparationActor_for_$jobTag"
    val jobPrepProps = JobPreparationActor.props(executionData, jobDescriptorKey, factory, dockerHashActor, initializationData, serviceRegistryActor, ioActor, backendSingletonActor)
    val jobPreparationActor = createJobPreparationActor(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! CallPreparation.Start
    goto(PreparingJob)
  }

  def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity): Try[ActorRef] = {
    val maybeFileHashingActorProps = factory.fileHashingActorProps map {
      _.apply(jobDescriptor, initializationData, serviceRegistryActor, ioActor)
    }

    maybeFileHashingActorProps match {
      case Some(fileHashingActorProps) =>
        val props = EngineJobHashingActor.props(
          self,
          jobDescriptor,
          initializationData,
          fileHashingActorProps,
          CallCacheReadingJobActor.props(callCacheReadActor),
          factory.runtimeAttributeDefinitions(initializationData), backendName, activity)
        val ejha = context.actorOf(props, s"ejha_for_$jobDescriptor")
        
        Success(ejha)
      case None => Failure(new IllegalStateException("Tried to initialize job hashing without a file hashing actor !"))
    }
  }

  def makeFetchCachedResultsActor(callCachingEntryId: CallCachingEntryId, taskOutputs: Seq[TaskOutput]): Unit = {
    context.actorOf(FetchCachedResultsActor.props(callCachingEntryId, self, new CallCache(SingletonServicesStore.databaseInterface)))
    ()
  }

  private def fetchCachedResults(taskOutputs: Seq[TaskOutput], callCachingEntryId: CallCachingEntryId, data: ResponsePendingData) = {
    makeFetchCachedResultsActor(callCachingEntryId, taskOutputs)
    goto(FetchingCachedOutputsFromDatabase) using data
  }

  private def makeBackendCopyCacheHit(wdlValueSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String], returnCode: Option[Int], data: ResponsePendingData, cacheResultId: CallCachingEntryId) = {
    factory.cacheHitCopyingActorProps match {
      case Some(propsMaker) =>
        val backendCacheHitCopyingActorProps = propsMaker(data.jobDescriptor, initializationData, serviceRegistryActor, ioActor)
        val cacheHitCopyActor = context.actorOf(backendCacheHitCopyingActorProps, buildCacheHitCopyingActorName(data.jobDescriptor, cacheResultId))
        cacheHitCopyActor ! CopyOutputsCommand(wdlValueSimpletons, jobDetritusFiles, returnCode)
        replyTo ! JobRunning(data.jobDescriptor.key, data.jobDescriptor.inputDeclarations, None)
        goto(BackendIsCopyingCachedOutputs)
      case None =>
        // This should be impossible with the FSM, but luckily, we CAN recover if some foolish future programmer makes this happen:
        val errorMessage = "Call caching copying should never have even been attempted with no copy actor props! (Programmer error!)"
        log.error(errorMessage)
        self ! JobFailedNonRetryableResponse(data.jobDescriptor.key, new RuntimeException(errorMessage), None)
        goto(BackendIsCopyingCachedOutputs)
    }
  }

  private def runJob(data: ResponsePendingData) = {
    val backendJobExecutionActor = context.actorOf(data.bjeaProps, buildJobExecutionActorName(data.jobDescriptor))
    val message = if (restarting) RecoverJobCommand else ExecuteJobCommand
    backendJobExecutionActor ! message
    replyTo ! JobRunning(data.jobDescriptor.key, data.jobDescriptor.inputDeclarations, Option(backendJobExecutionActor))
    goto(RunningJob) using data
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
        log.info("Trying to use another cache hit for job: {}", jobDescriptorKey)
        ejha ! NextHit
        goto(CheckingCallCache)
      case newData =>
        log.info("Could not find another cache hit, falling back to running job: {}", jobDescriptorKey)
        runJob(data)
    }
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowIdForLogging-BackendJobExecutionActor-$jobTag"
  }

  private def buildCacheHitCopyingActorName(jobDescriptor: BackendJobDescriptor, cacheResultId: CallCachingEntryId) = {
    s"$workflowIdForLogging-BackendCacheHitCopyingActor-$jobTag-${cacheResultId.id}"
  }

  private def invalidateCacheHitAndTransition(cacheId: CallCachingEntryId, data: ResponsePendingData, reason: Throwable) = {
    val invalidationRequired = effectiveCallCachingMode match {
      case CallCachingOff => throw new RuntimeException("Should not be calling invalidateCacheHit if call caching is off!") // Very unexpected. Fail out of this bad-state EJEA.
      case activity: CallCachingActivity => activity.options.invalidateBadCacheResults
    }
    if (invalidationRequired) {
      log.error(reason, "Failed copying cache results for job {}, invalidating cache entry.", jobDescriptorKey)
      invalidateCacheHit(cacheId)
      goto(InvalidatingCacheEntry)
    } else {
      handleCacheInvalidatedResponse(CallCacheInvalidationUnnecessary, data)
    }
  }

  protected def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = {
    val callCache = new CallCache(SingletonServicesStore.databaseInterface)
    context.actorOf(CallCacheInvalidateActor.props(callCache, cacheId), s"CallCacheInvalidateActor${cacheId.id}-$tag")
    ()
  }

  private def saveCacheResults(hashes: CallCacheHashes, data: SucceededResponseData) = {
    callCacheWriteActor ! SaveCallCacheHashes(workflowIdForLogging, hashes, data)
    val updatedData = data.copy(hashes = Option(Success(hashes)))
    goto(UpdatingCallCache) using updatedData
  }

  private def saveJobCompletionToJobStore(updatedData: ResponseData) = {
    updatedData.response match {
      case JobSucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: CallOutputs, _, executionEvents) =>
        saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case AbortedResponse(jobKey: BackendJobDescriptorKey) =>
        log.debug("{}: Won't save aborted job response to JobStore", jobTag)
        forwardAndStop(updatedData.response)
      case JobFailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case JobFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }
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

  private def addHashesAndStay(data: ResponsePendingData, hashes: CallCacheHashes): State = {
    val updatedData = data.copy(hashes = Option(Success(hashes)))
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
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState
  case object UpdatingCallCache extends EngineJobExecutionActorState
  case object UpdatingJobStore extends EngineJobExecutionActorState
  case object InvalidatingCacheEntry extends EngineJobExecutionActorState

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  val CacheMetadataKeyPrefix = CallMetadataKeys.CallCaching + MetadataKey.KeySeparator

  def props(replyTo: ActorRef,
            jobDescriptorKey: BackendJobDescriptorKey,
            executionData: WorkflowExecutionActorData,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            restarting: Boolean,
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef,
            jobStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            callCacheWriteActor: ActorRef,
            dockerHashActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonActor: Option[ActorRef],
            backendName: String,
            callCachingMode: CallCachingMode) = {
    Props(new EngineJobExecutionActor(
      replyTo = replyTo,
      jobDescriptorKey = jobDescriptorKey,
      executionData = executionData,
      factory = factory,
      initializationData = initializationData,
      restarting = restarting,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      jobStoreActor = jobStoreActor,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      dockerHashActor = dockerHashActor,
      jobTokenDispenserActor = jobTokenDispenserActor,
      backendSingletonActor = backendSingletonActor,
      backendName = backendName: String,
      callCachingMode = callCachingMode)).withDispatcher(EngineDispatcher)
  }

  private[execution] sealed trait EJEAData {
    override def toString = getClass.getSimpleName
  }

  private[execution] case object NoData extends EJEAData

  private[execution] case class ResponsePendingData(jobDescriptor: BackendJobDescriptor,
                                                    bjeaProps: Props,
                                                    hashes: Option[Try[CallCacheHashes]] = None,
                                                    ejha: Option[ActorRef] = None,
                                                    cacheHit: Option[CacheHit] = None) extends EJEAData {
    
    def withEJHA(ejha: ActorRef): EJEAData = this.copy(ejha = Option(ejha))


    def withSuccessResponse(success: JobSucceededResponse) = SucceededResponseData(success, hashes)

    def withResponse(response: BackendJobExecutionResponse) = response match {
      case success: JobSucceededResponse => SucceededResponseData(success, hashes)
      case failure => NotSucceededResponseData(failure, hashes)
    }

    def withCacheHit(cacheHit: CacheHit) = this.copy(cacheHit = Option(cacheHit))
  }

  private[execution] trait ResponseData extends EJEAData {
    def response: BackendJobExecutionResponse
    def hashes: Option[Try[CallCacheHashes]]
  }

  private[execution] case class SucceededResponseData(successResponse: JobSucceededResponse,
                                                      hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData {
    override def response = successResponse
  }

  private[execution] case class NotSucceededResponseData(response: BackendJobExecutionResponse,
                                                         hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData
}
