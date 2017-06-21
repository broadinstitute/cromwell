package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, RuntimeAttributeDefinition}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.core.callcaching._
import cromwell.core.logging.JobLogging
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.CallCachingHashEntry
import cromwell.engine.workflow.lifecycle.execution.CallMetadataHelper
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, FinalFileHashingResult, InitialHashingResult, NoFileHashesResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.services.metadata.CallMetadataKeys

/**
  * Coordinates the CallCacheHashingJobActor and the CallCacheReadingJobActor.
  * Over time will emit up to two messages back to its parent:
  *  * (if read enabled): Either a CacheHit(id) or CacheMiss message
  *  * (if write enabled): A CallCacheHashes(hashes) message
  */
class EngineJobHashingActor(receiver: ActorRef,
                            override val serviceRegistryActor: ActorRef,
                            jobDescriptor: BackendJobDescriptor,
                            initializationData: Option[BackendInitializationData],
                            fileHashingActorProps: Props,
                            callCacheReadingJobActorProps: Props,
                            runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                            backendName: String,
                            activity: CallCachingActivity,
                            callCachingEligible: CallCachingEligible) extends Actor with ActorLogging with JobLogging with CallMetadataHelper {

  override val jobTag = jobDescriptor.key.tag
  override val workflowId = jobDescriptor.workflowDescriptor.id
  override val workflowIdForCallMetadata: WorkflowId = workflowId

  private [callcaching] var initialHash: Option[InitialHashingResult] = None

  private [callcaching] val callCacheReadingJobActor = if (activity.readFromCache) {
    Option(context.actorOf(callCacheReadingJobActorProps))
  } else None

  override def preStart(): Unit = {
    context.actorOf(CallCacheHashingJobActor.props(
      jobDescriptor,
      callCacheReadingJobActor,
      initializationData,
      runtimeAttributeDefinitions,
      backendName,
      fileHashingActorProps,
      activity.writeToCache,
      callCachingEligible
    ))
    super.preStart()
  }

  override def receive = {
    case initialHashResult: InitialHashingResult => initialHash = Option(initialHashResult)
    case finalFileHashResult: FinalFileHashingResult => sendHashes(finalFileHashResult)
    case CacheMiss => receiver ! CacheMiss
    case hit: CacheHit => receiver ! hit
    case NextHit =>
      callCacheReadingJobActor match {
        case Some(readActor) =>
          readActor ! NextHit
        case None =>
          failAndStop(new IllegalStateException("Requested cache hit but there is no cache read actor"))
      }
    case hashingFailed: HashingFailedMessage =>
      failAndStop(hashingFailed.reason)
    case unexpected =>
      jobLogger.error(s"Received unexpected event $unexpected")
  }

  private def publishHashFailure(failure: Throwable) = {
    import cromwell.services.metadata.MetadataService._
    val failureAsEvents = throwableToMetadataEvents(metadataKeyForCall(jobDescriptor.key, CallMetadataKeys.CallCachingKeys.HashFailuresKey), failure)
    serviceRegistryActor ! PutMetadataAction(failureAsEvents)
  }

  private def failAndStop(reason: Throwable) = {
    publishHashFailure(reason)
    receiver ! HashError(reason)
    context stop self
  }

  private def sendHashes(finalFileHashingResult: FinalFileHashingResult) = {
    val fileHashes = finalFileHashingResult match {
      case CompleteFileHashingResult(fileHashResults, aggregatedFileHash) =>
        Option(FileHashes(fileHashResults, aggregatedFileHash))
      case NoFileHashesResult => None
    }

    initialHash match {
      case Some(initData) =>
        receiver ! CallCacheHashes(initData.initialHashes, initData.aggregatedBaseHash, fileHashes)
      case None =>
        failAndStop(new IllegalStateException("Received file hashes without initial hash."))
    }
  }
}

object EngineJobHashingActor {
  sealed trait EJHAState
  case object Running extends EJHAState
  case object WaitingForHashes extends EJHAState
  case object WaitingForJobSuccess extends EJHAState
  case object Done extends EJHAState

  object EJHAData {
    def empty = EJHAData(None)
  }

  case class EJHAData(initialHash: Option[InitialHashingResult])

  sealed trait EJHAResponse
  case object CacheMiss extends EJHAResponse
  case class CacheHit(cacheResultId: CallCachingEntryId) extends EJHAResponse
  case class HashError(reason: Throwable) extends EJHAResponse
  case class FileHashes(hashes: Set[HashResult], aggregatedHash: String)
  case class CallCacheHashes(initialHashes: Set[HashResult], aggregatedInitialHash: String, fileHashes: Option[FileHashes]) extends EJHAResponse {
    val hashes: Set[HashResult] = initialHashes ++ fileHashes.map(_.hashes).getOrElse(Set.empty)

    private lazy val hashPairs = hashes map {
      case HashResult(hashKey, HashValue(value)) => hashKey.key -> value
    }

    /**
      * Checks that the call caching join holds the same hash pairs as this instance
      */
    def equalsJoin(callCachingJoin: CallCachingJoin): Boolean = {

      val joinHashPairs = callCachingJoin.callCachingHashEntries.toSet[CallCachingHashEntry] map {
        case CallCachingHashEntry(key, value, _, _) => key -> value
      }

      (joinHashPairs, callCachingJoin.callCachingAggregationEntry) match {
        case (joinHashes, Some(joinAggregation))
          if joinAggregation.baseAggregation == aggregatedInitialHash &&
            joinAggregation.inputFilesAggregation == fileHashes.map(_.aggregatedHash) &&
            joinHashes.size == hashPairs.size &&
            joinHashes.diff(hashPairs).isEmpty =>
          true
        case _ => false
      }
    }
  }

  def props(receiver: ActorRef,
            serviceRegistryActor: ActorRef,
            jobDescriptor: BackendJobDescriptor,
            initializationData: Option[BackendInitializationData],
            fileHashingActorProps: Props,
            callCacheReadingJobActorProps: Props,
            runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
            backendName: String,
            activity: CallCachingActivity,
            callCachingEligible: CallCachingEligible): Props = Props(new EngineJobHashingActor(
    receiver = receiver,
    serviceRegistryActor = serviceRegistryActor,
    jobDescriptor = jobDescriptor,
    initializationData = initializationData,
    fileHashingActorProps = fileHashingActorProps,
    callCacheReadingJobActorProps = callCacheReadingJobActorProps,
    runtimeAttributeDefinitions = runtimeAttributeDefinitions,
    backendName = backendName,
    activity = activity,
    callCachingEligible = callCachingEligible)).withDispatcher(EngineDispatcher)
}
