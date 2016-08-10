package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.BackendSpecificHasherActor.{JobFileHashRequests, RuntimeAttributesHashesRequest, SingleFileHashRequest}
import wdl4s.values.WdlFile
import HashValue.StringMd5er

/**
  * Over time will emit up to two messages back to its parent:
  *  * (if read enabled): Either a CacheHit(id) or CacheMiss message
  *  * (if write enabled): A CallCacheHashes(hashes) message
  */
case class EngineJobHashingActor(jobDescriptor: BackendJobDescriptor,
                                 backendSpecificHasherActor: ActorRef,
                                 backendName: String,
                                 mode: CallCachingActivity) extends LoggingFSM[EJHAState, EJHAData] with ActorLogging {

  val receiver = context.parent

  initializeEJHA()

  when(DeterminingHitOrMiss) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) => lookupRelevantCacheResults(hashResultMessage.hashes)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) => checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResultIds(newCacheResults))
    case Event(PlaceholderHashKeyExpansion(placeholderHashKey, newHashKeys), _) => checkWhetherHitOrMissIsKnownThenTransition(stateData.replacePlaceholderHashKey(placeholderHashKey, newHashKeys))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
    case Event(PlaceholderHashKeyExpansion(placeholderHashKey, newHashKeys), _) => checkWhetherAllHashesAreKnownAndTransition(stateData.replacePlaceholderHashKey(placeholderHashKey, newHashKeys))
    case Event(CacheResultMatchesForHashes(_, _), _) => stay // Don't care; we already know the hit/miss status. Ignore this message
  }

  whenUnhandled {
    case Event(HashingFailedMessage(hashKey, reason), _) =>
      receiver ! HashError(new Exception(s"Unable to generate ${hashKey.key} hash. Caused by ${reason.getMessage}", reason))
      context.stop(self)
      stay
  }

  private def initializeEJHA() = {

    import cromwell.core.simpleton.WdlValueSimpleton._
    val inputSimpletons = jobDescriptor.inputs.simplify
    val (fileInputSimpletons, nonFileInputSimpletons) = inputSimpletons partition {
      case WdlValueSimpleton(_, f: WdlFile) => true
      case _ => false
    }

    val initialHashes = calculateInitialHashes(nonFileInputSimpletons, fileInputSimpletons)

    val fileContentHashesNeeded = if (mode.hashFileContents) {
      JobFileHashRequests(jobDescriptor.key, fileInputSimpletons map { case WdlValueSimpleton(name, x: WdlFile) => SingleFileHashRequest(HashKey(s"input: File(Content) $name"), x) })
    } else {
      JobFileHashRequests(jobDescriptor.key, List.empty)
    }
    val runtimeAttributeHashesNeededPlaceholder = HashKey("PLACEHOLDER: runtime attributes")

    val hashesNeeded: List[HashKey] = initialHashes.map(_.hashKey) ++ fileContentHashesNeeded.files.map(_.hashKey) :+ runtimeAttributeHashesNeededPlaceholder

    val initialState = if (mode.readFromCache) DeterminingHitOrMiss else GeneratingAllHashes
    val initialData = EJHAData(jobDescriptor, hashesNeeded)

    startWith(initialState, initialData)

    // Submit the set of initial hashes for checking against the DB:
    self ! EJHAInitialHashingResults(initialHashes)
    // Find the hashes for all input files:
    if (fileContentHashesNeeded.files.nonEmpty) backendSpecificHasherActor ! fileContentHashesNeeded
    backendSpecificHasherActor ! RuntimeAttributesHashesRequest(runtimeAttributeHashesNeededPlaceholder, jobDescriptor)
  }

  private def calculateInitialHashes(nonFileInputs: Iterable[WdlValueSimpleton], fileInputs: Iterable[WdlValueSimpleton]): List[HashResult] = {

    // TODO: Can we get instantiatedCommand in the JobPreparationActor?
    val commandTemplateHash = HashResult(HashKey("command template"), jobDescriptor.call.task.commandTemplate.toString.md5HashValue)
    val backendNameHash = HashResult(HashKey("backend name"), backendName.md5HashValue)

    val inputHashResults = nonFileInputs map {
      case WdlValueSimpleton(name, value) => HashResult(HashKey(s"input: ${value.wdlType.toWdlString} $name"),  value.toWdlString.md5HashValue)
    }

    val outputExpressionHashResults = jobDescriptor.call.task.outputs map { case output =>
      HashResult(HashKey(s"output expression: ${output.wdlType.toWdlString} ${output.name}"), output.requiredExpression.valueString.md5HashValue)
    }

    val inputFilePathHashes = if (mode.hashFilePaths) {
      fileInputs map { case WdlValueSimpleton(name, file: WdlFile) =>
        HashResult(HashKey(s"input: File(Path) $name"), file.value.toString.md5HashValue)
      }
    } else {
      List.empty[HashResult]
    }

    // Build these all together for the final set of initial hashes:
    List(commandTemplateHash, backendNameHash) ++ inputHashResults ++ outputExpressionHashResults ++ inputFilePathHashes
  }

  private def handleHashFailure(hashKey: HashKey) = {
    // TODO: do something better!
    self ! EJHAInitialHashingResults(List(HashResult(hashKey, HashValue("ERROR"))))
    stay
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
    val hitOrMissResponse = stateData.cacheHit match {
      case Some(cacheResultId) => CacheHit(cacheResultId)
      case None => CacheMiss
    }
    receiver ! hitOrMissResponse
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
}

object EngineJobHashingActor {

  def props(jobDescriptor: BackendJobDescriptor, fileHasherActor: ActorRef, backendName: String, activity: CallCachingActivity): Props = Props(new EngineJobHashingActor(jobDescriptor, fileHasherActor, backendName, activity))

  trait SuccessfulHashResultMessage {
    def hashes: Iterable[HashResult]
  }
  private[callcaching] case class EJHAInitialHashingResults(hashes: Iterable[HashResult]) extends SuccessfulHashResultMessage
  private[callcaching] case class HashingFailedMessage(key: HashKey, reason: Throwable)

  /**
    * Placeholder hash keys are used when we know the EJHA wants some hashes, but we don't know all the keys yet.
    *
    * E.g. We know we want runtime attributes hashes, but only the backend knows how many it will eventually create. So, we allow
    * the placeholder to be replaced via this expansion. And then we wait for a hash result for each of those new hash keys.
    */
  case class PlaceholderHashKeyExpansion(placeholderHashKey: HashKey, newHashKeysToAwait: Iterable[HashKey])

  sealed trait EJHAState
  case object DeterminingHitOrMiss extends EJHAState
  case object GeneratingAllHashes extends EJHAState

  sealed trait EJHAResponse
  case class CacheHit(cacheResultId: Int) extends EJHAResponse
  case object CacheMiss
  case class HashError(t: Throwable) {
    override def toString = s"HashError(${t.getMessage})"
  }
  case class CallCacheHashes(hashes: List[HashResult])
}

/**
  * Transient data for the EJHA.
  *
  * @param possibleCacheResults The set of cache results which have matched all currently tried hashes
  * @param cacheResultsIntersected The set of hash keys which have already been compared against the cache result table
  * @param hashesKnown The set of all hashes calculated so far (including initial hashes)
  * @param hashesNeeded The set of hashes we know that we need. Can be expanded via placeholder expansion
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

  def replacePlaceholderHashKey(placeholderHashKey: HashKey, newHashKeys: Iterable[HashKey]) = {
    val newHashesNeeded = this.hashesNeeded ++ newHashKeys filter {
      case `placeholderHashKey` => false // We don't want this!
      case _ => true // We do want all the others!
    }
    this.copy(hashesNeeded = newHashesNeeded)
  }

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
