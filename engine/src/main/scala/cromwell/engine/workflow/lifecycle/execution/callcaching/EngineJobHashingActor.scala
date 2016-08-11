package cromwell.engine.workflow.lifecycle.execution.callcaching

import java.nio.file.Path

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.database.sql.MetaInfoId
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.FileHasherActor.JobFileHashRequests

/**
  * Over time will emit up to two messages back to its parent:
  *  * (if read enabled): Either a CacheHit(id) or CacheMiss message
  *  * (if write enabled): A CallCacheHashes(hashes) message
  */
case class EngineJobHashingActor(jobDescriptor: BackendJobDescriptor,
                                 fileHasherActor: ActorRef,
                                 mode: CallCachingActivity) extends LoggingFSM[EJHAState, EJHAData] with ActorLogging {

  val receiver = context.parent

  initializeEJHA()

  when(DeterminingHitOrMiss) {
    case Event(hashResultMessage: HashResultMessage, _) =>
      // This is in DeterminingHitOrMiss, so use the new hash results to search for cache results.
      // Also update the state data with these new hash results.
      val newHashResults = hashResultMessage.hashes
      findCacheResults(newHashResults)
      stay using updateStateDataWithNewHashResults(newHashResults)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) =>
      checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResults(newCacheResults))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: HashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
    case Event(CacheResultMatchesForHashes(_, _), _) => stay // Don't care; we already know the hit/miss status. Ignore this message
  }

  whenUnhandled {
    case Event(failure: CacheResultLookupFailure, _) =>
      // Crash and let the supervisor deal with it.
      throw new RuntimeException("Failure looking up call cache results", failure.reason)
  }

  private def initializeEJHA() = {

    // TODO: Implement call caching hashing: #1230
    val initialHashes = Set(
      HashResult(HashKey("commandTemplate"), HashValue(jobDescriptor.call.task.commandTemplate.toString))
    )
    val fileHashesNeeded: Map[HashKey, Path] = Map.empty
    val hashesNeeded: Set[HashKey] = initialHashes.map(_.hashKey) ++ fileHashesNeeded.keys

    val initialState = if (mode.readFromCache) DeterminingHitOrMiss else GeneratingAllHashes
    val initialData = EJHAData(jobDescriptor, hashesNeeded)

    startWith(initialState, initialData)

    // Submit the set of initial hashes:
    self ! EJHAInitialHashingResults(initialHashes)
    // Find the hashes for all input files:
    fileHasherActor ! allInputJobFileHashRequests
  }

  private def checkWhetherHitOrMissIsKnownThenTransition(newData: EJHAData) = {
    if (newData.isDefinitelyCacheHitOrMiss) {
      respondWithHitOrMissThenTransition(newData)
    }
    else {
      stay() using newData // Stay in DeterminingHitOrMiss
    }
  }

  private def respondWithHitOrMissThenTransition(newData: EJHAData) = {
    newData.cacheHit match {
      case Some(cacheResultId) => receiver ! CacheHit(cacheResultId)
      case None => receiver ! CacheMiss
    }
    if (!mode.writeToCache) {
      context.stop(self)
      stay
    } else {
      checkWhetherAllHashesAreKnownAndTransition(newData)
    }
  }

  private def checkWhetherAllHashesAreKnownAndTransition(newData: EJHAData) = {
    if (newData.allHashesKnown) {
      receiver ! CallCacheHashes(newData.hashesKnown)
      context.stop(self)
    }
    goto(GeneratingAllHashes) using newData
  }

  /**
    * Needs to convert a hash result into the set of CachedResults which are consistent with it
    */
  private def findCacheResults(hashResults: Set[HashResult]) = {
    val hashes = CallCacheHashes(hashResults)
    context.actorOf(CallCacheReadActor.props(hashes), s"CallCacheReadActor-${jobDescriptor.descriptor.id}-${jobDescriptor.key.tag}")
  }

  def updateStateDataWithNewHashResults(hashResults: Set[HashResult]): EJHAData = {
    if (mode.writeToCache) stateData.withNewKnownHashes(hashResults) else stateData
  }

  // TODO: Implement me.
  // Called at actor start-time. Should convert all job input files into file hash requests
  private def allInputJobFileHashRequests: JobFileHashRequests = JobFileHashRequests(jobDescriptor.key, List.empty)
}

object EngineJobHashingActor {

  def props(jobDescriptor: BackendJobDescriptor, fileHasherActor: ActorRef, activity: CallCachingActivity): Props = Props(new EngineJobHashingActor(jobDescriptor, fileHasherActor, activity))

  trait HashResultMessage {
    def hashes: Set[HashResult]
  }
  private[callcaching] case class EJHAInitialHashingResults(hashes: Set[HashResult]) extends HashResultMessage
  private[callcaching] case object CheckWhetherAllHashesAreKnown

  sealed trait EJHAState
  case object DeterminingHitOrMiss extends EJHAState
  case object GeneratingAllHashes extends EJHAState

  sealed trait EJHAResponse
  case class CacheHit(cacheResultId: MetaInfoId) extends EJHAResponse
  case object CacheMiss
  case class CallCacheHashes(hashes: Set[HashResult])
}

/**
  * Transient data for the EJHA.
  *
  * @param possibleCacheResults The set of cache results which have matched all currently tried hashes
  * @param keysCheckedAgainstCache The set of hash keys which have already been compared against the cache result table
  * @param hashesKnown The set of all hashes calculated so far (including initial hashes)
  * @param hashesNeeded Not transient but oh-so-useful for everything else.
  */
private[callcaching] case class EJHAData(possibleCacheResults: Set[MetaInfoId],
                                         keysCheckedAgainstCache: Set[HashKey],
                                         hashesKnown: Set[HashResult],
                                         hashesNeeded: Set[HashKey]) {

  // Manipulators
  def intersectCacheResults(newCacheResults: CacheResultMatchesForHashes) = {
    val newIds = newCacheResults.cacheResultIds
    val intersectedIds = if (keysCheckedAgainstCache.nonEmpty) possibleCacheResults.intersect(newIds) else newIds
    this.copy(
      possibleCacheResults = intersectedIds,
      keysCheckedAgainstCache = keysCheckedAgainstCache ++ newCacheResults.hashResults.map(_.hashKey))
  }
  def withNewKnownHashes(hashResults: Set[HashResult]) = this.copy(hashesKnown = hashesKnown ++ hashResults) // NB doesn't update hashesNeeded which should be considered immutable

  // Queries
  def allHashesKnown = hashesNeeded.size == hashesKnown.size
  def allCacheResultsIntersected = keysCheckedAgainstCache.size == hashesNeeded.size
  def cacheHit = if (allCacheResultsIntersected && possibleCacheResults.nonEmpty) possibleCacheResults.headOption else None
  def isDefinitelyCacheHit = cacheHit.isDefined
  def isDefinitelyCacheMiss = keysCheckedAgainstCache.nonEmpty && possibleCacheResults.isEmpty
  def isDefinitelyCacheHitOrMiss = isDefinitelyCacheHit || isDefinitelyCacheMiss
}

object EJHAData {
  def apply(jobDescriptor: BackendJobDescriptor, hashesNeeded: Set[HashKey]): EJHAData = EJHAData(Set.empty, Set.empty, Set.empty, hashesNeeded)
}
