package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.database.CromwellDatabase
import cromwell.database.sql.MetaInfoId
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, CallCacheHashes, HashError}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.{Pending => _, _}
import cromwell.core.ExecutionIndex.IndexEnhancedIndex
import scala.util.{Failure, Success, Try}

class EngineJobExecutionActor(jobKey: BackendJobDescriptorKey,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              serviceRegistryActor: ActorRef,
                              jobStoreActor: ActorRef,
                              backendName: String,
                              callCachingMode: CallCachingMode) extends LoggingFSM[EngineJobExecutionActorState, EJEAData] with WorkflowLogging {

  override val workflowId = executionData.workflowDescriptor.id

  val jobTag = s"${workflowId.shortString}:${jobKey.call.fullyQualifiedName}:${jobKey.index.fromIndex}:${jobKey.attempt}"
  val tag = s"EJEA_$jobTag"

  // There's no need to check for a cache hit again if we got preempted
  // NB: this can also change (e.g. if we have a HashError we just force this to CallCachingOff)
  private var effectiveCallCachingMode = if (jobKey.attempt > 1) { callCachingMode.withoutRead } else { callCachingMode }

  val effectiveCallCachingKey = "Effective call caching mode"

  log.info(s"$tag: $effectiveCallCachingKey: ${effectiveCallCachingMode.getClass.getSimpleName}")
  writeToMetadata(Map(effectiveCallCachingKey -> effectiveCallCachingMode.getClass.getSimpleName))

  startWith(Pending, NoData)

  // When Pending, the FSM always has NoData
  when(Pending) {
    case Event(Execute, NoData) =>
      if (restarting) {
        val jobStoreKey = jobKey.toJobStoreKey(workflowId)
        jobStoreActor ! QueryJobCompletion(jobStoreKey)
        goto(CheckingJobStore)
      } else {
        prepareJob()
      }
  }

  // When CheckingJobStore, the FSM always has NoData
  when(CheckingJobStore) {
    case Event(JobNotComplete, NoData) =>
      prepareJob()
    case Event(JobComplete(jobResult), NoData) =>
      jobResult match {
        case JobResultSuccess(returnCode, jobOutputs) =>
          context.parent ! SucceededResponse(jobKey, returnCode, jobOutputs)
          context stop self
          stay()
        case JobResultFailure(returnCode, reason, false) =>
          context.parent ! FailedNonRetryableResponse(jobKey, reason, returnCode)
          context stop self
          stay()
        case JobResultFailure(returnCode, reason, true) =>
          context.parent ! FailedRetryableResponse(jobKey, reason, returnCode)
          context stop self
          stay()
      }
    case Event(f: JobStoreReadFailure, NoData) =>
      log.error(f.reason, "{}: Error reading from JobStore", tag)
      // Escalate
      throw new RuntimeException(f.reason)
  }

  // When PreparingJob, the FSM always has NoData
  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, bjeaProps), NoData) =>
      effectiveCallCachingMode match {
        case activity: CallCachingActivity if activity.readFromCache =>
          initializeJobHashing(jobDescriptor, activity)
          goto(CheckingCallCache) using JobDescriptorData(jobDescriptor, bjeaProps)
        case activity: CallCachingActivity =>
          initializeJobHashing(jobDescriptor, activity)
          runJob(jobDescriptor, bjeaProps)
        case CallCachingOff => runJob(jobDescriptor, bjeaProps)
      }
    case Event(response: BackendJobPreparationFailed, NoData) =>
      context.parent forward response
      context stop self
      stay()
  }

  // When CheckingCallCache, the FSM always has EJEAJobDescriptorData
  private val callCachingReadResultMetadataKey = "Call caching read result"
  when(CheckingCallCache) {
    case Event(HashError(t), JobDescriptorData(jobDescriptor, bjeaProps)) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Hashing Error: ${t.getMessage}"))
      recordHashError(t)
      runJob(jobDescriptor, bjeaProps)
    case Event(CacheMiss, JobDescriptorData(jobDescriptor, bjeaProps)) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> "Cache Miss"))
      log.info(s"Cache miss for job ${jobDescriptor.key.call.fullyQualifiedName}, index ${jobDescriptor.key.index}")
      runJob(jobDescriptor, bjeaProps)
    case Event(CacheHit(cacheResultId), JobDescriptorData(jobDescriptor, _)) =>
      writeToMetadata(Map(callCachingReadResultMetadataKey -> s"Cache Hit (from result ID $cacheResultId)"))
      log.info(s"Cache hit for job ${jobDescriptor.key.call.fullyQualifiedName}, index ${jobDescriptor.key.index}! Copying cache result $cacheResultId")
      lookupCachedResult(jobDescriptor, cacheResultId)
  }

  // When RunningJob, the FSM always has EJEAPartialCompletionData (which might be None, None)
  when(RunningJob) {
    case Event(hashes: CallCacheHashes, data @ PartialCompletionDataWithSucceededResponse(response: SucceededResponse)) =>
      saveCacheResults(CacheWriteOnCompletionData(response, hashes))
    case Event(hashes: CallCacheHashes, data @ EmptyPartialCompletionData) =>
      stay using PartialCompletionDataWithHashes(Success(hashes))

    case Event(HashError(t), data @ PartialCompletionDataWithSucceededResponse(response)) =>
      recordHashError(t)
      saveJobCompletionToJobStore(response)
    case Event(HashError(t), data @ EmptyPartialCompletionData) =>
      recordHashError(t)
      stay using PartialCompletionDataWithHashes(Failure(t))

    case Event(response: SucceededResponse, PartialCompletionDataWithHashes(Success(hashes))) if effectiveCallCachingMode.writeToCache =>
      saveCacheResults(CacheWriteOnCompletionData(response, hashes))
    case Event(response: SucceededResponse, data @ EmptyPartialCompletionData) if effectiveCallCachingMode.writeToCache =>
      stay using PartialCompletionDataWithSucceededResponse(response)
    case Event(response: BackendJobExecutionResponse, _) =>
      saveJobCompletionToJobStore(response)
  }

  // When WritingToCallCache, the FSM always has EJEASuccessfulCompletionDataWithHashes
  when(UpdatingCallCache) {
    case Event(CallCacheWriteSuccess, data @ CacheWriteOnCompletionData(response, _)) =>
      saveJobCompletionToJobStore(response)
    case Event(CallCacheWriteFailure(reason), data @ CacheWriteOnCompletionData(response, _)) =>
      context.parent ! FailedNonRetryableResponse(jobKey, reason, response.returnCode)
      context stop self
      stay()
  }

  // When UpdatingJobStore, the FSM always has EJEACompletionData
  when(UpdatingJobStore) {
    case Event(JobStoreWriteSuccess(_), CacheWriteOffCompletionData(response)) =>
      context.parent forward response
      context stop self
      stay()
    case Event(JobStoreWriteFailure(t), CacheWriteOffCompletionData(_)) =>
      context.parent ! FailedNonRetryableResponse(jobKey, new Exception(s"JobStore write failure: ${t.getMessage}", t), None)
      context.stop(self)
      stay()
  }

  onTransition {
    case fromState -> toState =>
      log.debug("Transitioning from {}({}) to {}({})", fromState, stateData, toState, nextStateData)
  }

  whenUnhandled {
    case Event(msg, _) =>
      log.error(s"Bad message to EngineJobExecutionActor in state $stateName(with data $stateData): $msg")
      stay
  }

  private def recordHashError(reason: Throwable) = {
    log.error("{}: hash error: {}. Disabling call caching for this job", tag, reason.getMessage)
    effectiveCallCachingMode = CallCachingOff
  }

  def prepareJob() = {
    val jobPreparationActorName = s"BackendPreparationActor_for_$jobTag"
    val jobPrepProps = JobPreparationActor.props(executionData, jobKey, factory, initializationData, serviceRegistryActor)
    val jobPreparationActor = context.actorOf(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! JobPreparationActor.Start
    goto(PreparingJob)
  }

  def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity) = {
    val fileHasherActor = context.actorOf(BackendSpecificHasherActor.props(activity),  s"FileHasherActor_for_$jobTag")
    context.actorOf(EngineJobHashingActor.props(jobDescriptor, fileHasherActor, backendName, activity))
  }

  def lookupCachedResult(jobDescriptor: BackendJobDescriptor, cacheResultId: MetaInfoId) = {
    // TODO: Start up a backend job copying actor (if possible, otherwise just runJob). That should send back the BackendJobExecutionResponse
    self ! FailedNonRetryableResponse(jobKey, new Exception("Call cache result copying not implemented!"), None)
    // While the cache result is looked up, we wait for the response just like we were waiting for a Job to complete:
    goto(RunningJob) using EmptyPartialCompletionData
  }

  def runJob(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) = {
    val backendJobExecutionActor = context.actorOf(bjeaProps, buildJobExecutionActorName(jobDescriptor))
    val message = if (restarting) RecoverJobCommand else ExecuteJobCommand
    backendJobExecutionActor ! message
    context.parent ! JobRunning(jobDescriptor, backendJobExecutionActor)
    goto(RunningJob) using EmptyPartialCompletionData
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }

  private def saveCacheResults(completionData: CacheWriteOnCompletionData) = {
    val callCache = new CallCache(CromwellDatabase.databaseInterface)
    context.actorOf(CallCacheWriteActor.props(callCache, workflowId, completionData.hashes, completionData.jobResult), s"CallCacheWriteActor-$tag")
    goto(UpdatingCallCache) using completionData
  }

  private def saveJobCompletionToJobStore(response: BackendJobExecutionResponse) = {
    response match {
      case SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs) => saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case AbortedResponse(jobKey: BackendJobDescriptorKey) => log.debug("Won't save 'aborted' job response to JobStore")
      case FailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case FailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }
    goto(UpdatingJobStore) using CacheWriteOffCompletionData(response)
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
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoputter
    serviceRegistryActor.putMetadata(workflowId, Option(jobKey), keyValues)
  }
}

object EngineJobExecutionActor {
  /** States */
  sealed trait EngineJobExecutionActorState
  case object Pending extends EngineJobExecutionActorState
  case object CheckingJobStore extends EngineJobExecutionActorState
  case object CheckingCallCache extends EngineJobExecutionActorState
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState
  case object UpdatingCallCache extends EngineJobExecutionActorState
  case object UpdatingJobStore extends EngineJobExecutionActorState

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  final case class JobRunning(jobDescriptor: BackendJobDescriptor, backendJobExecutionActor: ActorRef)

  def props(jobDescriptorKey: BackendJobDescriptorKey,
            executionData: WorkflowExecutionActorData,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            restarting: Boolean,
            serviceRegistryActor: ActorRef,
            jobStoreActor: ActorRef,
            backendName: String,
            callCachingMode: CallCachingMode) = {
    Props(new EngineJobExecutionActor(jobDescriptorKey,
      executionData,
      factory,
      initializationData,
      restarting,
      serviceRegistryActor,
      jobStoreActor,
      backendName: String,
      callCachingMode)).withDispatcher(EngineDispatcher)
  }
}

private[execution] sealed trait EJEAData { override def toString = getClass.getSimpleName }

private[execution] case object NoData extends EJEAData
private[execution] case class JobDescriptorData(jobDescriptor: BackendJobDescriptor, bjeaActorProps: Props) extends EJEAData

private[execution] case object EmptyPartialCompletionData extends EJEAData
private[execution] case class PartialCompletionDataWithSucceededResponse(response: SucceededResponse) extends EJEAData
private[execution] case class PartialCompletionDataWithHashes(hashes: Try[CallCacheHashes]) extends EJEAData

private[execution] case class CacheWriteOffCompletionData(jobResult: BackendJobExecutionResponse) extends EJEAData
private[execution] case class CacheWriteOnCompletionData(jobResult: SucceededResponse, hashes: CallCacheHashes) extends EJEAData
