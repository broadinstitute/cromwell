package cromwell.engine.workflow.lifecycle.execution.callcaching

import java.nio.file.Path

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobDescriptor
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
    case Event(hashResultMessage: HashResultMessage, _) => lookupRelevantCacheResults(hashResultMessage.hashes)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) => checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResultIds(newCacheResults))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: HashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
    case Event(CacheResultMatchesForHashes(_, _), _) => stay // Don't care; we already know the hit/miss status. Ignore this message
  }

  private def initializeEJHA() = {

    // TODO: Implement call caching hashing: #1230
    val initialHashes = List(
      HashResult(HashKey("commandTemplate"), HashValue(jobDescriptor.call.task.commandTemplate.toString))
    )
    val fileHashesNeeded: Map[HashKey, Path] = Map.empty
    val hashesNeeded: List[HashKey] = initialHashes.map(_.hashKey) ++ fileHashesNeeded.keys

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
    stateData.cacheHit match {
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
  private def lookupRelevantCacheResults(hashResults: Iterable[HashResult]) = {
    val newData = if (mode.writeToCache) stateData.withNewKnownHashes(hashResults) else stateData
    self ! CacheResultMatchesForHashes(hashResults, Set.empty) // TODO: Probably want to convert this into a call to another actor which would respond back with the results.
    stay using newData
  }

  // TODO: Implement me.
  // Called at actor start-time. Should convert all job input files into file hash requests
  private def allInputJobFileHashRequests: JobFileHashRequests = JobFileHashRequests(jobDescriptor.key, List.empty)
}

object EngineJobHashingActor {

  def props(jobDescriptor: BackendJobDescriptor, fileHasherActor: ActorRef, activity: CallCachingActivity): Props = Props(new EngineJobHashingActor(jobDescriptor, fileHasherActor, activity))

  trait HashResultMessage {
    def hashes: Iterable[HashResult]
  }
  private[callcaching] case class EJHAInitialHashingResults(hashes: Iterable[HashResult]) extends HashResultMessage
  private[callcaching] case object CheckWhetherAllHashesAreKnown

  sealed trait EJHAState
  case object DeterminingHitOrMiss extends EJHAState
  case object GeneratingAllHashes extends EJHAState

  sealed trait EJHAResponse
  case class CacheHit(cacheResultId: Int) extends EJHAResponse
  case object CacheMiss
  case class CallCacheHashes(hashes: List[HashResult])
}

/**
  * Transient data for the EJHA.
  *
  * @param possibleCacheResults The set of cache results which have matched all currently tried hashes
  * @param cacheResultsIntersected The set of hash keys which have already been compared against the cache result table
  * @param hashesKnown The set of all hashes calculated so far (including initial hashes)
  * @param hashesNeeded Not transient but oh-so-useful for everything else.
  */
private[callcaching] case class EJHAData(possibleCacheResults: Set[Int],
                                         cacheResultsIntersected: List[HashKey],
                                         hashesKnown: List[HashResult],
                                         hashesNeeded: Iterable[HashKey]) {

  // Manipulators
  def intersectCacheResultIds(newCacheResults: CacheResultMatchesForHashes) = {
    val newIds = newCacheResults.cacheResultIds
    val intersectedIds = if (cacheResultsIntersected.nonEmpty) possibleCacheResults.intersect(newIds) else newIds
    this.copy(
      possibleCacheResults = intersectedIds,
      cacheResultsIntersected = cacheResultsIntersected ++ newCacheResults.hashResults.map(_.hashKey))
  }
  def withNewKnownHashes(hashResults: Iterable[HashResult]) = this.copy(hashesKnown = hashesKnown ++ hashResults) // NB doesn't update hashesNeeded which should be considered immutable

  // Queries
  def allHashesKnown = hashesNeeded.size == hashesKnown.size
  def allCacheResultsIntersected = cacheResultsIntersected.size == hashesNeeded.size
  def cacheHit = if (allCacheResultsIntersected && possibleCacheResults.nonEmpty) Option(possibleCacheResults.head) else None
  def isDefinitelyCacheHit = cacheHit.isDefined
  def isDefinitelyCacheMiss = cacheResultsIntersected.nonEmpty && possibleCacheResults.isEmpty
  def isDefinitelyCacheHitOrMiss = isDefinitelyCacheHit || isDefinitelyCacheMiss
}

object EJHAData {
  def apply(jobDescriptor: BackendJobDescriptor, hashesNeeded: Iterable[HashKey]): EJHAData = EJHAData(Set.empty, List.empty, List.empty, hashesNeeded)
}
