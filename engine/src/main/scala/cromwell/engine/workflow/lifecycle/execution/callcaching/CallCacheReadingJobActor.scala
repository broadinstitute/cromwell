package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching.HashingFailedMessage
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCachePathPrefixes
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, HashError}
import cromwell.services.CallCaching.CallCachingEntryId

/**
  * Receives hashes from the CallCacheHashingJobActor and makes requests to the database to determine whether or not there might be a hit
  * for this job.
  * 
  * First receives the initial hashes, and asks the database if there is at least one entry with the same aggregated initial hash.
  * If not, it's a CacheMiss.
  * If yes, ask the CallCacheHashingJobActor for the next batch of file hashes.
  * Every time a new batch of hashes is received, check against the database if at least one entry matches all those hashes.
  * Keep asking for new batches until either one returns no matching entry, in which case it's a CacheMiss, or until it receives
  * the last batch along with the aggregated file hash.
  * In the latter case, asks the database for the first entry matching both the initial and aggregated file hash (if any).
  * Sends the response to its parent.
  * In case of a CacheHit, stays alive in case using the hit fails and it needs to fetch the next one. Otherwise just dies.
  */
class CallCacheReadingJobActor(callCacheReadActor: ActorRef, prefixesHint: Option[CallCachePathPrefixes]) extends LoggingFSM[CallCacheReadingJobActorState, CCRJAData] {
  
  startWith(WaitingForInitialHash, CCRJANoData)
  
  when(WaitingForInitialHash) {
    case Event(InitialHashingResult(_, aggregatedBaseHash, hints), CCRJANoData) =>
      callCacheReadActor ! HasMatchingInitialHashLookup(aggregatedBaseHash, hints)
      goto(WaitingForHashCheck) using CCRJAWithData(sender(), aggregatedBaseHash, fileHash = None, seenCaches = Set.empty)
  }
  
  when(WaitingForHashCheck) {
    case Event(HasMatchingEntries, CCRJAWithData(hashingActor, _, _, _)) =>
      hashingActor ! NextBatchOfFileHashesRequest
      goto(WaitingForFileHashes)
    case Event(NoMatchingEntries, _) =>
      cacheMiss
  }
  
  when(WaitingForFileHashes) {
    case Event(CompleteFileHashingResult(_, aggregatedFileHash), data: CCRJAWithData) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(data.initialHash, aggregatedFileHash), data.seenCaches, prefixesHint)
      goto(WaitingForCacheHitOrMiss) using data.withFileHash(aggregatedFileHash)
    case Event(NoFileHashesResult, data: CCRJAWithData) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(data.initialHash, None), data.seenCaches, prefixesHint)
      goto(WaitingForCacheHitOrMiss)
  }
  
  when(WaitingForCacheHitOrMiss) {
    case Event(CacheLookupNextHit(hit), data: CCRJAWithData) =>
      context.parent ! CacheHit(hit)
      stay() using data.withSeenCache(hit)
    case Event(CacheLookupNoHit, _) =>
      cacheMiss
    case Event(NextHit, CCRJAWithData(_, aggregatedInitialHash, aggregatedFileHash, seenCaches)) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, aggregatedFileHash), seenCaches, prefixesHint)
      stay()
  }

  whenUnhandled {
    case Event(_: HashingFailedMessage, _) =>
      // No need to send to the parent since it also receives file hash updates
      cacheMiss
    case Event(CacheResultLookupFailure(failure), _) =>
      context.parent ! HashError(failure)
      cacheMiss
  }

  def cacheMiss = {
    context.parent ! CacheMiss
    context stop self
    stay()
  }
}

object CallCacheReadingJobActor {

  def props(callCacheReadActor: ActorRef, prefixesHint: Option[CallCachePathPrefixes]) = {
    Props(new CallCacheReadingJobActor(callCacheReadActor, prefixesHint)).withDispatcher(EngineDispatcher)
  }

  sealed trait CallCacheReadingJobActorState
  case object WaitingForInitialHash extends CallCacheReadingJobActorState
  case object WaitingForHashCheck extends CallCacheReadingJobActorState
  case object WaitingForFileHashes extends CallCacheReadingJobActorState
  case object WaitingForCacheHitOrMiss extends CallCacheReadingJobActorState

  sealed trait CCRJAData
  case object CCRJANoData extends CCRJAData
  case class CCRJAWithData(hashingActor: ActorRef, initialHash: String, fileHash: Option[String], seenCaches: Set[CallCachingEntryId]) extends CCRJAData {
    def withSeenCache(id: CallCachingEntryId): CCRJAWithData = this.copy(seenCaches = seenCaches + id)
    def withFileHash(aggregatedFileHash: String): CCRJAWithData = this.copy(fileHash = Option(aggregatedFileHash))
  }

  case object NextHit
}
