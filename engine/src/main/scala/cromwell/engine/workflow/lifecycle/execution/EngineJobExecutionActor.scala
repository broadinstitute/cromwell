package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.database.CromwellDatabase
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, CallCacheHashes}
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.{Pending => _, _}

class EngineJobExecutionActor(jobKey: BackendJobDescriptorKey,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              serviceRegistryActor: ActorRef,
                              jobStoreActor: ActorRef,
                              callCachingMode: CallCachingMode) extends LoggingFSM[EngineJobExecutionActorState, EJEAData] with WorkflowLogging {

  override val workflowId = executionData.workflowDescriptor.id

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
      log.error(f.reason, "Error reading from JobStore for " + jobKey)
      // Escalate
      throw new RuntimeException(f.reason)
  }

  // When PreparingJob, the FSM always has NoData
  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, bjeaProps), NoData) =>
      callCachingMode.activity match {
        case Some(activity) if activity.readFromCache =>
          initializeJobHashing(jobDescriptor, activity)
          goto(CheckingCallCache) using EJEAJobDescriptorData(Option(jobDescriptor), Option(bjeaProps))
        case Some(activity) =>
          initializeJobHashing(jobDescriptor, activity)
          runJob(jobDescriptor, bjeaProps)
        case None => runJob(jobDescriptor, bjeaProps)
      }
    case Event(response: BackendJobPreparationFailed, NoData) =>
      context.parent forward response
      context stop self
      stay()
  }

  // When CheckingCallCache, the FSM always has EJEAJobDescriptorData
  when(CheckingCallCache) {
    case Event(CacheMiss, EJEAJobDescriptorData(Some(jobDescriptor), Some(bjeaProps))) =>
      log.info(s"Cache miss for job ${jobDescriptor.key.call.fullyQualifiedName}, index ${jobDescriptor.key.index}")
      runJob(jobDescriptor, bjeaProps)
    case Event(CacheHit(cacheResultId), EJEAJobDescriptorData(Some(jobDescriptor), _)) =>
      log.info(s"Cache hit for job ${jobDescriptor.key.call.fullyQualifiedName}, index ${jobDescriptor.key.index}! Copying cache result $cacheResultId")
      lookupCachedResult(jobDescriptor, cacheResultId)
  }

  // When RunningJob, the FSM always has EJEAPartialCompletionData (which might be None, None)
  when(RunningJob) {
    case Event(response: SucceededResponse, EJEAPartialCompletionData(None, Some(hashes))) if callCachingMode.writeToCache =>
      saveCacheResults(EJEASuccessfulCompletionDataWithHashes(response, hashes))
    case Event(response: SucceededResponse, data @ EJEAPartialCompletionData(None, None)) if callCachingMode.writeToCache =>
      stay using data.copy(jobResult = Option(response))
    case Event(hashes: CallCacheHashes, data @ EJEAPartialCompletionData(Some(response: SucceededResponse), None)) =>
      saveCacheResults(EJEASuccessfulCompletionDataWithHashes(response, hashes))
    case Event(hashes: CallCacheHashes, data @ EJEAPartialCompletionData(None, None)) =>
      stay using data.copy(hashes = Option(hashes))
    case Event(response: BackendJobExecutionResponse, data: EJEAPartialCompletionData) =>
      saveJobCompletionToJobStore(response)
  }

  // When WritingToCallCache, the FSM always has EJEASuccessfulCompletionDataWithHashes
  when(UpdatingCallCache) {
    case Event(CallCacheWriteSuccess, data @ EJEASuccessfulCompletionDataWithHashes(response, _)) =>
      saveJobCompletionToJobStore(response)
    case Event(CallCacheWriteFailure(reason), data @ EJEASuccessfulCompletionDataWithHashes(response, _)) =>
      context.parent ! FailedNonRetryableResponse(jobKey, reason, response.returnCode)
      context stop self
      stay()
  }

  // When UpdatingJobStore, the FSM always has EJEACompletionData
  when(UpdatingJobStore) {
    case Event(JobStoreWriteSuccess(_), EJEACompletionData(response)) =>
      context.parent forward response
      context stop self
      stay()
    case Event(JobStoreWriteFailure(t), EJEACompletionData(_)) =>
      context.parent ! FailedNonRetryableResponse(jobKey, new Exception(s"JobStore write failure: ${t.getMessage}", t), None)
      context.stop(self)
      stay()
  }

  onTransition {
    case fromState -> toState =>
      log.info(s"Transitioning from $fromState($stateData) to $toState($nextStateData)")
  }

  whenUnhandled {
    case Event(msg, _) =>
      log.error(s"Bad message to EngineJobExecutionActor in state $stateName(with data $stateData): $msg")
      stay
  }

  def prepareJob() = {
    val jobPreparationActorName = s"$workflowId-BackendPreparationActor-${jobKey.tag}"
    val jobPrepProps = JobPreparationActor.props(executionData, jobKey, factory, initializationData, serviceRegistryActor)
    val jobPreparationActor = context.actorOf(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! JobPreparationActor.Start
    goto(PreparingJob)
  }

  def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity) = {
    val fileHasherActor = context.actorOf(FileHasherActor.props)
    context.actorOf(EngineJobHashingActor.props(jobDescriptor, fileHasherActor, activity))
  }

  def lookupCachedResult(jobDescriptor: BackendJobDescriptor, cacheResultId: Int) = {
    // TODO: Start up a backend job copying actor (if possible, otherwise just runJob). That should send back the BackendJobExecutionResponse
    self ! FailedNonRetryableResponse(jobKey, new Exception("Call cache writing incomplete! Turn it off!!"), None)
    // While the cache result is looked up, we wait for the resoonse just like we were waiting for a Job to complete:
    goto(RunningJob) using EJEAPartialCompletionData(None, None)
  }

  def runJob(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) = {
    val backendJobExecutionActor = context.actorOf(bjeaProps, buildJobExecutionActorName(jobDescriptor))
    val message = if (restarting) RecoverJobCommand else ExecuteJobCommand
    backendJobExecutionActor ! message
    context.parent ! JobRunning(jobDescriptor, backendJobExecutionActor)
    goto(RunningJob) using EJEAPartialCompletionData(None, None)
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }

  private def saveCacheResults(completionData: EJEASuccessfulCompletionDataWithHashes) = {
    val callCache = new CallCache(CromwellDatabase.databaseInterface)
    context.actorOf(CallCacheWriteActor.props(callCache, workflowId, completionData.hashes, completionData.jobResult))
    goto(UpdatingCallCache) using completionData
  }

  private def saveJobCompletionToJobStore(response: BackendJobExecutionResponse) = {
    response match {
      case SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs) => saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
      case AbortedResponse(jobKey: BackendJobDescriptorKey) => log.debug("Won't save 'aborted' job response to JobStore")
      case FailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
      case FailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
    }
    goto(UpdatingJobStore) using EJEACompletionData(response)
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
            callCachingMode: CallCachingMode) = {
    Props(new EngineJobExecutionActor(jobDescriptorKey,
      executionData,
      factory,
      initializationData,
      restarting,
      serviceRegistryActor,
      jobStoreActor,
      callCachingMode)).withDispatcher(EngineDispatcher)
  }
}

sealed trait EJEAData
case object NoData extends EJEAData
case class EJEAJobDescriptorData(jobDescriptor: Option[BackendJobDescriptor], bjeaActorProps: Option[Props]) extends EJEAData {
  override def toString = s"EJEAJobDescriptorData(backendJobDescriptor: ${jobDescriptor.isDefined}, bjeaActorProps: ${bjeaActorProps.isDefined})"
}

case class EJEAPartialCompletionData(jobResult: Option[BackendJobExecutionResponse], hashes: Option[CallCacheHashes]) extends EJEAData {
  override def toString = s"EJEAPartialCompletionData(jobResult: ${jobResult.isDefined}, hashes: ${hashes.isDefined})"
}

// Le sigh. Can't do case-to-case inheritance so let's jump through these extractor-pattern hoops...
class EJEACompletionData(val jobResult: BackendJobExecutionResponse) extends EJEAData {
  override def toString = s"EJEACompletionData(${jobResult.getClass.getSimpleName})"
}
object EJEACompletionData { def apply(jobResult: BackendJobExecutionResponse) = new EJEACompletionData(jobResult); def unapply(completionData: EJEACompletionData) = Some(completionData.jobResult) }
case class EJEASuccessfulCompletionDataWithHashes(override val jobResult: SucceededResponse, hashes: CallCacheHashes) extends EJEACompletionData(jobResult = jobResult) {
  override def toString = s"EJEACompletionDataWithHashes(${jobResult.getClass.getSimpleName}, ${hashes.hashes.size} hashes)"
}
