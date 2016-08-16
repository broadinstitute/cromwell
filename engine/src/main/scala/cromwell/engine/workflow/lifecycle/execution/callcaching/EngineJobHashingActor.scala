package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.database.sql.MetaInfoId
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
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) =>
      // This is in DeterminingHitOrMiss, so use the new hash results to search for cache results.
      // Also update the state data with these new hash results.
      val newHashResults = hashResultMessage.hashes
      findCacheResults(newHashResults)
      stay using updateStateDataWithNewHashResults(newHashResults)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) =>
      checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResults(newCacheResults))
    case Event(RuntimeAttributesHashKeyPlaceholderExpansion(newHashKeys), _) => checkWhetherHitOrMissIsKnownThenTransition(stateData.replacePlaceholderHashKey(RuntimeAttributeHashKeyPlaceholder, newHashKeys))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
    case Event(RuntimeAttributesHashKeyPlaceholderExpansion(newHashKeys), _) => checkWhetherAllHashesAreKnownAndTransition(stateData.replacePlaceholderHashKey(RuntimeAttributeHashKeyPlaceholder, newHashKeys))
    case Event(CacheResultMatchesForHashes(_, _), _) => stay // Don't care; we already know the hit/miss status. Ignore this message
  }

  whenUnhandled {
    case Event(failure: CacheResultLookupFailure, _) =>
      // Crash and let the supervisor deal with it.
      throw new RuntimeException("Failure looking up call cache results", failure.reason)
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

    val fileContentHashesNeeded = mode.fileHashingType match {
      case HashFileContents => JobFileHashRequests(jobDescriptor.key, fileInputSimpletons map { case WdlValueSimpleton(name, x: WdlFile) => SingleFileHashRequest(HashKey(s"input: File(Content) $name"), x) })
      case HashFilePath => JobFileHashRequests(jobDescriptor.key, List.empty)
    }

    val hashesNeeded: Set[HashKey] = initialHashes.map(_.hashKey) ++ fileContentHashesNeeded.files.map(_.hashKey) ++ Set(RuntimeAttributeHashKeyPlaceholder)

    val initialState = if (mode.readFromCache) DeterminingHitOrMiss else GeneratingAllHashes
    val initialData = EJHAData(jobDescriptor, hashesNeeded)

    startWith(initialState, initialData)

    // Submit the set of initial hashes for checking against the DB:
    self ! EJHAInitialHashingResults(initialHashes)
    // Find the hashes for all input files:
    if (fileContentHashesNeeded.files.nonEmpty) backendSpecificHasherActor ! fileContentHashesNeeded
    backendSpecificHasherActor ! RuntimeAttributesHashesRequest(jobDescriptor)
  }

  private def calculateInitialHashes(nonFileInputs: Iterable[WdlValueSimpleton], fileInputs: Iterable[WdlValueSimpleton]): Set[HashResult] = {

    // TODO: Can we get instantiatedCommand in the JobPreparationActor?
    val commandTemplateHash = HashResult(HashKey("command template"), jobDescriptor.call.task.commandTemplate.toString.md5HashValue)
    val backendNameHash = HashResult(HashKey("backend name"), backendName.md5HashValue)

    val inputHashResults = nonFileInputs map {
      case WdlValueSimpleton(name, value) => HashResult(HashKey(s"input: ${value.wdlType.toWdlString} $name"),  value.toWdlString.md5HashValue)
    }

    val outputExpressionHashResults = jobDescriptor.call.task.outputs map { case output =>
      HashResult(HashKey(s"output expression: ${output.wdlType.toWdlString} ${output.name}"), output.requiredExpression.valueString.md5HashValue)
    }

    val inputFilePathHashes = mode.fileHashingType match {
      case HashFilePath =>
        fileInputs map {
          case WdlValueSimpleton(name, file: WdlFile) => HashResult(HashKey(s"input: File(Path) $name"), file.value.toString.md5HashValue)
        }
      case HashFileContents => List.empty[HashResult]
    }

    // Build these all together for the final set of initial hashes:
    Set(commandTemplateHash, backendNameHash) ++ inputHashResults ++ outputExpressionHashResults ++ inputFilePathHashes
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
    val hitOrMissResponse = newData.cacheHit match {
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
  private var lookupIndex = 0
  private def findCacheResults(hashResults: Set[HashResult]) = if ( hashResults.nonEmpty) {
      val hashes = CallCacheHashes(hashResults)
      val actorName = s"CallCacheReadActor-${jobDescriptor.descriptor.id}-${jobDescriptor.key.tag}-lookup_$lookupIndex"
      lookupIndex += 1
      context.actorOf(CallCacheReadActor.props(hashes), actorName)
  } else ()

  def updateStateDataWithNewHashResults(hashResults: Set[HashResult]): EJHAData = {
    if (mode.writeToCache) stateData.withNewKnownHashes(hashResults) else stateData
  }
}

object EngineJobHashingActor {

  def props(jobDescriptor: BackendJobDescriptor, fileHasherActor: ActorRef, backendName: String, activity: CallCachingActivity): Props = Props(new EngineJobHashingActor(jobDescriptor, fileHasherActor, backendName, activity))

  trait SuccessfulHashResultMessage {
    def hashes: Set[HashResult]
  }
  private[callcaching] case class EJHAInitialHashingResults(hashes: Set[HashResult]) extends SuccessfulHashResultMessage
  private[callcaching] case object CheckWhetherAllHashesAreKnown
  private[callcaching] case class HashingFailedMessage(key: HashKey, reason: Throwable)

  /**
    * Placeholder hash keys are used when we know the EJHA wants some hashes, but we don't know all the keys yet.
    *
    * E.g. We know we want runtime attributes hashes, but only the backend knows how many it will eventually create. So, we allow
    * the placeholder to be replaced via this expansion. And then we wait for a hash result for each of those new hash keys.
    */
  case class RuntimeAttributesHashKeyPlaceholderExpansion(newHashKeysToAwait: Iterable[HashKey])

  sealed trait EJHAState
  case object DeterminingHitOrMiss extends EJHAState
  case object GeneratingAllHashes extends EJHAState

  sealed trait EJHAResponse
  case class CacheHit(cacheResultId: MetaInfoId) extends EJHAResponse
  case object CacheMiss
  case class HashError(t: Throwable) {
    override def toString = s"HashError(${t.getMessage})"
  }
  case class CallCacheHashes(hashes: Set[HashResult])
}

/**
  * Transient data for the EJHA.
  *
  * @param possibleCacheResults The set of cache results which have matched all currently tried hashes
  * @param keysCheckedAgainstCache The set of hash keys which have already been compared against the cache result table
  * @param hashesKnown The set of all hashes calculated so far (including initial hashes)
  * @param hashesNeeded The set of hashes we know that we need. Can be expanded via placeholder expansion
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
  def withNewKnownHashes(hashResults: Set[HashResult]) = this.copy(hashesKnown = hashesKnown ++ hashResults)

  def replacePlaceholderHashKey(placeholderHashKey: HashKey, newHashKeys: Iterable[HashKey]) = {
    val newHashesNeeded = this.hashesNeeded ++ newHashKeys filter {
      case `placeholderHashKey` => false // We d  on't want this!
      case _ => true // We do want all the others!
    }
    this.copy(hashesNeeded = newHashesNeeded)
  }

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
