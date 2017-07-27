package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.callcaching.HashingFailedMessage
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult, PartialFileHashingResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, HashError}
import cromwell.core.Dispatcher.EngineDispatcher

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
class CallCacheReadingJobActor(callCacheReadActor: ActorRef) extends LoggingFSM[CallCacheReadingJobActorState, CCRJAData] {
  
  startWith(WaitingForInitialHash, CCRJANoData)
  
  when(WaitingForInitialHash) {
    case Event(InitialHashingResult(_, aggregatedBaseHash), CCRJANoData) =>
      callCacheReadActor ! HasMatchingInitialHashLookup(aggregatedBaseHash)
      goto(WaitingForHashCheck) using CCRJAWithData(sender(), aggregatedBaseHash, None, 1)
  }
  
  when(WaitingForHashCheck) {
    case Event(HasMatchingEntries, CCRJAWithData(hashingActor, _, _, _)) =>
      hashingActor ! NextBatchOfFileHashesRequest
      goto(WaitingForFileHashes)
    case Event(NoMatchingEntries, _) =>
      cacheMiss
  }
  
  when(WaitingForFileHashes) {
    case Event(PartialFileHashingResult(hashes), _) =>
      callCacheReadActor ! HasMatchingInputFilesHashLookup(hashes)
      goto(WaitingForHashCheck)
    case Event(CompleteFileHashingResult(_, aggregatedFileHash), data: CCRJAWithData) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(data.initialHash, aggregatedFileHash), data.currentHitNumber)
      goto(WaitingForCacheHitOrMiss) using data.withFileHash(aggregatedFileHash)
    case Event(NoFileHashesResult, data: CCRJAWithData) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(data.initialHash, None), data.currentHitNumber)
      goto(WaitingForCacheHitOrMiss)
  }
  
  when(WaitingForCacheHitOrMiss) {
    case Event(CacheLookupNextHit(hit), data: CCRJAWithData) =>
      context.parent ! CacheHit(hit)
      stay() using data.increment
    case Event(CacheLookupNoHit, _) =>
      cacheMiss
    case Event(NextHit, CCRJAWithData(_, aggregatedInitialHash, aggregatedFileHash, currentHitNumber)) =>
      callCacheReadActor ! CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, aggregatedFileHash), currentHitNumber)
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
  
  def props(callCacheReadActor: ActorRef) = {
    Props(new CallCacheReadingJobActor(callCacheReadActor)).withDispatcher(EngineDispatcher)
  }
  
  sealed trait CallCacheReadingJobActorState
  case object WaitingForInitialHash extends CallCacheReadingJobActorState
  case object WaitingForHashCheck extends CallCacheReadingJobActorState
  case object WaitingForFileHashes extends CallCacheReadingJobActorState
  case object WaitingForCacheHitOrMiss extends CallCacheReadingJobActorState

  sealed trait CCRJAData
  case object CCRJANoData extends CCRJAData
  case class CCRJAWithData(hashingActor: ActorRef, initialHash: String, fileHash: Option[String], currentHitNumber: Int) extends CCRJAData {
    def increment = this.copy(currentHitNumber = currentHitNumber + 1)
    def withFileHash(aggregatedFileHash: String) = this.copy(fileHash = Option(aggregatedFileHash))
  }

  case object NextHit
}
