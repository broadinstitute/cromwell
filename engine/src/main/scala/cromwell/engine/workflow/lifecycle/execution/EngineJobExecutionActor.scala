package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import akka.actor.{ActorRef, ActorRefFactory, LoggingFSM, Props}
import akka.routing.RoundRobinPool
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, CallCacheHashes, HashError}
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDenied, JobExecutionTokenDispensed, JobExecutionTokenRequest, JobExecutionTokenReturn}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.{Pending => _, _}
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataJobKey, MetadataKey, MetadataValue}
import wdl4s.TaskOutput

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success, Try}

class EngineJobExecutionActor(replyTo: ActorRef,
                              jobDescriptorKey: BackendJobDescriptorKey,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              serviceRegistryActor: ActorRef,
                              jobStoreActor: ActorRef,
                              callCacheReadActor: ActorRef,
                              jobTokenDispenserActor: ActorRef,
                              backendSingletonActor: Option[ActorRef],
                              backendName: String,
                              callCachingMode: CallCachingMode) extends LoggingFSM[EngineJobExecutionActorState, EJEAData] with WorkflowLogging {

  override val workflowId = executionData.workflowDescriptor.id

  val jobTag = s"${workflowId.shortString}:${jobDescriptorKey.call.fullyQualifiedName}:${jobDescriptorKey.index.fromIndex}:${jobDescriptorKey.attempt}"
  val tag = s"EJEA_$jobTag"

  // There's no need to check for a cache hit again if we got preempted, or if there's no result copying actor defined
  // NB: this can also change (e.g. if we have a HashError we just force this to CallCachingOff)
  private var effectiveCallCachingMode = if (factory.cacheHitCopyingActorProps.isEmpty || jobDescriptorKey.attempt > 1) callCachingMode.withoutRead else callCachingMode

  // For tests:
  private[execution] def checkEffectiveCallCachingMode = effectiveCallCachingMode

  private[execution] var executionToken: Option[JobExecutionToken] = None

  private val effectiveCallCachingKey = "Effective call caching mode"

  implicit val ec: ExecutionContext = context.dispatcher

  log.debug(s"$tag: $effectiveCallCachingKey: $effectiveCallCachingMode")
  writeCallCachingModeToMetadata()

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
      if (restarting) {
        val jobStoreKey = jobDescriptorKey.toJobStoreKey(workflowId)
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
        case JobResultSuccess(returnCode, jobOutputs) => SucceededResponse(jobDescriptorKey, returnCode, jobOutputs, None, Seq.empty)
        case JobResultFailure(returnCode, reason, false) => FailedNonRetryableResponse(jobDescriptorKey, reason, returnCode)
        case JobResultFailure(returnCode, reason, true) => FailedRetryableResponse(jobDescriptorKey, reason, returnCode)
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
        case activity: CallCachingActivity if activity.readFromCache =>
          initializeJobHashing(jobDescriptor, activity)
          goto(CheckingCallCache) using updatedData
        case activity: CallCachingActivity =>
          initializeJobHashing(jobDescriptor, activity)
          runJob(updatedData)
        case CallCachingOff => runJob(updatedData)
      }
    case Event(response: BackendJobPreparationFailed, NoData) =>
      forwardAndStop(response)
  }

  private val callCachingReadResultMetadataKey = "Call caching read result"
  when(CheckingCallCache) {
    case Event(CacheMiss, data: ResponsePendingData) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> "Cache Miss"))
      log.debug("Cache miss for job {}", jobTag)
      runJob(data)
    case Event(hit @ CacheHit(cacheResultId), data: ResponsePendingData) =>
      fetchCachedResults(data, jobDescriptorKey.call.task.outputs, hit)
    case Event(HashError(t), data: ResponsePendingData) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Hashing Error: ${t.getMessage}"))
      disableCallCaching(t)
      runJob(data)
  }

  when(FetchingCachedOutputsFromDatabase) {
    case Event(CachedOutputLookupSucceeded(wdlValueSimpletons, jobDetritus, returnCode, cacheResultId, cacheHitDetails), data: ResponsePendingData) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Cache Hit: $cacheHitDetails"))
      log.debug("Cache hit for {}! Fetching cached result {}", jobTag, cacheResultId)
      makeBackendCopyCacheHit(cacheResultId, wdlValueSimpletons, jobDetritus, returnCode, data)
    case Event(CachedOutputLookupFailed(metaInfoId, error), data: ResponsePendingData) =>
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
    case Event(response: SucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)))) =>
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: SucceededResponse, data @ ResponsePendingData(_, _, None)) if effectiveCallCachingMode.writeToCache =>
      // Wait for the CallCacheHashes
      stay using data.withSuccessResponse(response)
    case Event(response: SucceededResponse, data: ResponsePendingData) => // bad hashes or cache write off
      saveJobCompletionToJobStore(data.withSuccessResponse(response))
    case Event(response: BackendJobExecutionResponse, data: ResponsePendingData) =>
      // This matches all response types other than `SucceededResponse`.
      response match {
        case f: BackendJobFailedResponse =>log.error("{}: Failed copying cache results, falling back to running job: {}", jobDescriptorKey, f.throwable)
        case _ => //
      }
      runJob(data)

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

  when(RunningJob) {
    case Event(hashes: CallCacheHashes, data: SucceededResponseData) =>
      saveCacheResults(hashes, data)
    case Event(hashes: CallCacheHashes, data: ResponsePendingData) =>
      addHashesAndStay(data, hashes)

    case Event(HashError(t), data: SucceededResponseData) =>
      disableCallCaching(t)
      saveJobCompletionToJobStore(data.copy(hashes = Option(Failure(t))))
    case Event(HashError(t), data: ResponsePendingData) =>
      disableCallCaching(t)
      stay using data.copy(hashes = Option(Failure(t)))

    case Event(response: SucceededResponse, data @ ResponsePendingData(_, _, Some(Success(hashes)))) if effectiveCallCachingMode.writeToCache =>
      eventList ++= response.executionEvents
      saveCacheResults(hashes, data.withSuccessResponse(response))
    case Event(response: SucceededResponse, data @ ResponsePendingData(_, _, None)) if effectiveCallCachingMode.writeToCache =>
      log.debug(s"Got job result for {}, awaiting hashes", jobTag)
      stay using data.withSuccessResponse(response)
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
      respondAndStop(FailedNonRetryableResponse(jobDescriptorKey, new Exception(s"JobStore write failure: ${t.getMessage}", t), None))
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

  private def requestExecutionToken(): Unit = {
    jobTokenDispenserActor ! JobExecutionTokenRequest(factory.jobExecutionTokenType)
  }

  private def returnExecutionToken(): Unit = {
    executionToken foreach { jobTokenDispenserActor ! JobExecutionTokenReturn(_) }
  }

  private def forwardAndStop(response: Any): State = {
    replyTo forward response
    returnExecutionToken()
    tellEventMetadata()
    context stop self
    stay()
  }

  private def respondAndStop(response: Any): State = {
    replyTo ! response
    returnExecutionToken()
    tellEventMetadata()
    context stop self
    stay()
  }

  private def disableCallCaching(reason: Throwable) = {
    log.error(reason, "{}: Hash error, disabling call caching for this job.", jobTag)
    effectiveCallCachingMode = CallCachingOff
    writeCallCachingModeToMetadata()
  }

  private def disableCacheWrite(reason: Throwable) = {
    log.error("{}: Disabling cache writing for this job.", jobTag)
    if (effectiveCallCachingMode.writeToCache) {
      effectiveCallCachingMode = effectiveCallCachingMode.withoutWrite
      writeCallCachingModeToMetadata()
    }
  }

  def writeCallCachingModeToMetadata(): Unit = {
    writeToMetadata(Map(effectiveCallCachingKey -> effectiveCallCachingMode.toString))
  }

  def createJobPreparationActor(jobPrepProps: Props, name: String): ActorRef = context.actorOf(jobPrepProps, name)
  def prepareJob() = {
    val jobPreparationActorName = s"BackendPreparationActor_for_$jobTag"
    val jobPrepProps = JobPreparationActor.props(executionData, jobDescriptorKey, factory, initializationData, serviceRegistryActor, backendSingletonActor)
    val jobPreparationActor = createJobPreparationActor(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! JobPreparationActor.Start
    goto(PreparingJob)
  }

  def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity): Unit = {
    val props = EngineJobHashingActor.props(
      self,
      jobDescriptor,
      initializationData,
      // Use context.system instead of context as the factory. Otherwise when we die, so will the child actors.
      factoryFileHashingRouter(backendName, factory, context.system),
      callCacheReadActor,
      factory.runtimeAttributeDefinitions(initializationData), backendName, activity)
    context.actorOf(props, s"ejha_for_$jobDescriptor")
    ()
  }

  def makeFetchCachedResultsActor(cacheHit: CacheHit, taskOutputs: Seq[TaskOutput]): Unit = {
    context.actorOf(FetchCachedResultsActor.props(cacheHit, self, new CallCache(SingletonServicesStore.databaseInterface)))
    ()
  }

  private def fetchCachedResults(data: ResponsePendingData, taskOutputs: Seq[TaskOutput], cacheHit: CacheHit) = {
    makeFetchCachedResultsActor(cacheHit, taskOutputs)
    goto(FetchingCachedOutputsFromDatabase)
  }

  private def makeBackendCopyCacheHit(cacheHit: CacheHit, wdlValueSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String], returnCode: Option[Int], data: ResponsePendingData) = {
    factory.cacheHitCopyingActorProps match {
      case Some(propsMaker) =>
        val backendCacheHitCopyingActorProps = propsMaker(data.jobDescriptor, initializationData, serviceRegistryActor)
        val cacheHitCopyActor = context.actorOf(backendCacheHitCopyingActorProps, buildCacheHitCopyingActorName(data.jobDescriptor))
        cacheHitCopyActor ! CopyOutputsCommand(wdlValueSimpletons, jobDetritusFiles, returnCode)
        replyTo ! JobRunning(data.jobDescriptor, None)
        goto(BackendIsCopyingCachedOutputs) using data
      case None =>
        // This should be impossible with the FSM, but luckily, we CAN recover if some foolish future programmer makes this happen:
        val errorMessage = "Call caching copying should never have even been attempted with no copy actor props! (Programmer error!)"
        log.error(errorMessage)
        self ! FailedNonRetryableResponse(data.jobDescriptor.key, new RuntimeException(errorMessage), None)
        goto(BackendIsCopyingCachedOutputs)
    }
  }

  private def runJob(data: ResponsePendingData) = {
    val backendJobExecutionActor = context.actorOf(data.bjeaProps, buildJobExecutionActorName(data.jobDescriptor))
    val message = if (restarting) RecoverJobCommand else ExecuteJobCommand
    backendJobExecutionActor ! message
    replyTo ! JobRunning(data.jobDescriptor, Option(backendJobExecutionActor))
    goto(RunningJob) using data
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-$jobTag"
  }

  private def buildCacheHitCopyingActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendCacheHitCopyingActor-$jobTag"
  }

  protected def createSaveCacheResultsActor(hashes: CallCacheHashes, success: SucceededResponse): Unit = {
    val callCache = new CallCache(SingletonServicesStore.databaseInterface)
    context.actorOf(CallCacheWriteActor.props(callCache, workflowId, hashes, success), s"CallCacheWriteActor-$tag")
    ()
  }

  private def saveCacheResults(hashes: CallCacheHashes, data: SucceededResponseData) = {
    createSaveCacheResultsActor(hashes, data.successResponse)
    val updatedData = data.copy(hashes = Option(Success(hashes)))
    goto(UpdatingCallCache) using updatedData
  }

  private def saveJobCompletionToJobStore(updatedData: ResponseData) = {
    updatedData.response match {
      case SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs, _, _) => saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case AbortedResponse(jobKey: BackendJobDescriptorKey) =>
        log.debug("{}: Won't save aborted job response to JobStore", jobTag)
        forwardAndStop(updatedData.response)
      case FailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case FailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }
    goto(UpdatingJobStore) using updatedData
  }

  private def saveSuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], outputs: JobOutputs) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowId)
    val jobStoreResult = JobResultSuccess(returnCode, outputs)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }

  private def saveUnsuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], reason: Throwable, retryable: Boolean) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowId)
    val jobStoreResult = JobResultFailure(returnCode, reason, retryable)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }

  private def writeToMetadata(keyValues: Map[String, String]) = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(workflowId, Option(jobDescriptorKey), keyValues)
  }

  private def addHashesAndStay(data: ResponsePendingData, hashes: CallCacheHashes): State = {
    val updatedData = data.copy(hashes = Option(Success(hashes)))
    stay using updatedData
  }

  /**
    * Fire and forget events to the metadata service
    */
  private def tellEventMetadata(): Unit = {
    eventList.headOption foreach { firstEvent =>
      // The final event is only used as the book-end for the final pairing so the name is never actually used...
      val offset = firstEvent.offsetDateTime.getOffset
      val now = OffsetDateTime.now.withOffsetSameInstant(offset)
      val lastEvent = ExecutionEvent("!!Bring Back the Monarchy!!", now)
      val tailedEventList = eventList :+ lastEvent
      val events = tailedEventList.sliding(2).zipWithIndex flatMap {
        case (Seq(eventCurrent, eventNext), index) =>
          val eventKey = s"executionEvents[$index]"
          List(
            metadataEvent(s"$eventKey:description", eventCurrent.name),
            metadataEvent(s"$eventKey:startTime", eventCurrent.offsetDateTime),
            metadataEvent(s"$eventKey:endTime", eventNext.offsetDateTime)
          )
      }

      serviceRegistryActor ! PutMetadataAction(events.toIterable)
    }
  }

  private def metadataEvent(key: String, value: Any) = {
    val metadataValue = MetadataValue(value)
    MetadataEvent(metadataKey(key), metadataValue)
  }

  private lazy val metadataJobKey = {
    MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)
  }
  private def metadataKey(key: String) = MetadataKey(workflowId, Option(metadataJobKey), key)
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

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  final case class JobRunning(jobDescriptor: BackendJobDescriptor, backendJobExecutionActor: Option[ActorRef])

  def props(replyTo: ActorRef,
            jobDescriptorKey: BackendJobDescriptorKey,
            executionData: WorkflowExecutionActorData,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            restarting: Boolean,
            serviceRegistryActor: ActorRef,
            jobStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
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
      jobStoreActor = jobStoreActor,
      callCacheReadActor = callCacheReadActor,
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
                                                    hashes: Option[Try[CallCacheHashes]] = None) extends EJEAData {

    def withSuccessResponse(success: SucceededResponse) = SucceededResponseData(success, hashes)

    def withResponse(response: BackendJobExecutionResponse) = response match {
      case success: SucceededResponse => SucceededResponseData(success, hashes)
      case failure => NotSucceededResponseData(failure, hashes)
    }
  }

  private[execution] trait ResponseData extends EJEAData {
    def response: BackendJobExecutionResponse
    def hashes: Option[Try[CallCacheHashes]]
  }

  private[execution] case class SucceededResponseData(successResponse: SucceededResponse,
                                                      hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData {
    override def response = successResponse
  }

  private[execution] case class NotSucceededResponseData(response: BackendJobExecutionResponse,
                                                         hashes: Option[Try[CallCacheHashes]] = None) extends ResponseData

  /**
    * Deliberately a singleton (well, a singleton router), so we can globally rate limit hash lookups per backend.
    *
    * More refinement may appear via #1377.
    */
  private var factoryFileHashingRouters = Map[BackendLifecycleActorFactory, ActorRef]()

  /**
    * Returns a RoundRobinPool of actors based on the backend factory.
    *
    * @param backendName                  Name of the backend.
    * @param backendLifecycleActorFactory A backend factory.
    * @param actorRefFactory              An actor factory.
    * @return a RoundRobinPool of actors based on backend factory.
    */
  private def factoryFileHashingRouter(backendName: String,
                                       backendLifecycleActorFactory: BackendLifecycleActorFactory,
                                       actorRefFactory: ActorRefFactory): ActorRef = {
    synchronized {
      val (originalOrUpdated, result) = getOrElseUpdated(
        factoryFileHashingRouters, backendLifecycleActorFactory, {
          val numberOfInstances = backendLifecycleActorFactory.fileHashingActorCount
          val props = backendLifecycleActorFactory.fileHashingActorProps
          actorRefFactory.actorOf(RoundRobinPool(numberOfInstances).props(props), s"FileHashingActor-$backendName")
        }
      )
      factoryFileHashingRouters = originalOrUpdated
      result
    }
  }

  /**
    * Immutable version of mutable.Map.getOrElseUpdate based on:
    * http://stackoverflow.com/questions/4385976/idiomatic-get-or-else-update-for-immutable-map#answer-5840119
    *
    * If given key is already in this map, returns associated value in the copy of the Map.
    *
    * Otherwise, computes value from given expression `op`, stores with key
    * in map and returns that value in a copy of the Map.
    *
    * @param map the immutable map
    * @param key the key to test
    * @param op  the computation yielding the value to associate with `key`, if
    *            `key` is previously unbound.
    * @tparam K type of the key
    * @tparam V type of the value
    * @return the value associated with key (either previously or as a result
    *         of executing the method).
    */
  def getOrElseUpdated[K, V](map: Map[K, V], key: K, op: => V): (Map[K, V], V) = {
    map.get(key) match {
      case Some(value) => (map, value)
      case None =>
        val value = op
        (map.updated(key, value), value)
    }
  }
}
