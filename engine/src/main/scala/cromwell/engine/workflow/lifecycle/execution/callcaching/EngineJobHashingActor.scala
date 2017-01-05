package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, RuntimeAttributeDefinition}
import cromwell.core.callcaching._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.{CacheLookupRequest, CacheResultLookupFailure, CacheResultMatchesForHashes}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import wdl4s.values.WdlFile

/**
  * Over time will emit up to two messages back to its parent:
  *  * (if read enabled): Either a CacheHit(id) or CacheMiss message
  *  * (if write enabled): A CallCacheHashes(hashes) message
  */
case class EngineJobHashingActor(receiver: ActorRef,
                                 jobDescriptor: BackendJobDescriptor,
                                 initializationData: Option[BackendInitializationData],
                                 fileHashingActor: ActorRef,
                                 callCacheReadActor: ActorRef,
                                 runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                                 backendName: String,
                                 activity: CallCachingActivity) extends LoggingFSM[EJHAState, EJHAData] with ActorLogging {

  initializeEJHA()

  when(DeterminingHitOrMiss) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) =>
      // This is in DeterminingHitOrMiss, so use the new hash results to search for cache results.
      // Also update the state data with these new hash results.
      val newHashResults = hashResultMessage.hashes
      findCacheResults(newHashResults)
      updateStateDataWithNewHashResultsAndTransition(newHashResults)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) =>
      checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResults(newCacheResults))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
    case Event(CacheResultMatchesForHashes(_, _), _) => stay // Don't care; we already know the hit/miss status. Ignore this message
  }

  whenUnhandled {
    case Event(CacheResultLookupFailure(reason), _) =>
      receiver ! HashError(new Exception(s"Failure looking up call cache results: ${reason.getMessage}"))
      context.stop(self)
      stay
    case Event(HashingFailedMessage(hashKey, reason), _) =>
      receiver ! HashError(new Exception(s"Unable to generate ${hashKey.key} hash. Caused by ${reason.getMessage}", reason))
      context.stop(self)
      stay
    case Event(other, _) =>
      log.error(s"Bad message in $stateName with $stateData: $other")
      stay
  }

  onTransition {
    case fromState -> toState =>
      log.debug("Transitioning from {}({}) to {}({})", fromState, stateData, toState, nextStateData)
  }

  private def initializeEJHA() = {

    import cromwell.core.simpleton.WdlValueSimpleton._

    val inputSimpletons = jobDescriptor.fullyQualifiedInputs.simplify
    val (fileInputSimpletons, nonFileInputSimpletons) = inputSimpletons partition {
      case WdlValueSimpleton(_, f: WdlFile) => true
      case _ => false
    }

    val initialHashes = calculateInitialHashes(nonFileInputSimpletons, fileInputSimpletons)

    val fileContentHashesNeeded = fileInputSimpletons collect {
      case WdlValueSimpleton(name, x: WdlFile) => SingleFileHashRequest(jobDescriptor.key, HashKey(s"input: File $name"), x, initializationData)
    }

    val hashesNeeded: Set[HashKey] = initialHashes.map(_.hashKey) ++ fileContentHashesNeeded.map(_.hashKey)

    val initialState = if (activity.readFromCache) DeterminingHitOrMiss else GeneratingAllHashes
    val initialData = EJHAData(hashesNeeded, activity)

    startWith(initialState, initialData)

    // Submit the set of initial hashes for checking against the DB:
    self ! EJHAInitialHashingResults(initialHashes)
    // Find the hashes for all input files:
    fileContentHashesNeeded.foreach(fileHashingActor ! _)
  }

  private def calculateInitialHashes(nonFileInputs: Iterable[WdlValueSimpleton], fileInputs: Iterable[WdlValueSimpleton]): Set[HashResult] = {

    val commandTemplateHash = HashResult(HashKey("command template"), jobDescriptor.call.task.commandTemplateString.md5HashValue)
    val backendNameHash = HashResult(HashKey("backend name"), backendName.md5HashValue)
    val inputCountHash = HashResult(HashKey("input count"), (nonFileInputs.size + fileInputs.size).toString.md5HashValue)
    val outputCountHash = HashResult(HashKey("output count"), jobDescriptor.call.task.outputs.size.toString.md5HashValue)

    val runtimeAttributeHashes = runtimeAttributeDefinitions map { definition => jobDescriptor.runtimeAttributes.get(definition.name) match {
      case Some(wdlValue) => HashResult(HashKey("runtime attribute: " + definition.name, definition.usedInCallCaching), wdlValue.valueString.md5HashValue)
      case None => HashResult(HashKey("runtime attribute: " + definition.name, definition.usedInCallCaching), UnspecifiedRuntimeAttributeHashValue)
    }}

    val inputHashResults = nonFileInputs map {
      case WdlValueSimpleton(name, value) => HashResult(HashKey(s"input: ${value.wdlType.toWdlString} $name"),  value.toWdlString.md5HashValue)
    }

    val outputExpressionHashResults = jobDescriptor.call.task.outputs map { output =>
      HashResult(HashKey(s"output expression: ${output.wdlType.toWdlString} ${output.unqualifiedName}"), output.requiredExpression.valueString.md5HashValue)
    }

    // Build these all together for the final set of initial hashes:
    Set(commandTemplateHash, backendNameHash, inputCountHash, outputCountHash) ++ runtimeAttributeHashes ++ inputHashResults ++ outputExpressionHashResults
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
    import cats.data.NonEmptyList
    val hitOrMissResponse: EJHAResponse = newData.cacheHits map { _.toList } flatMap NonEmptyList.fromList map CacheHit.apply getOrElse CacheMiss

    receiver ! hitOrMissResponse
    if (!activity.writeToCache) {
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
    val filtered = hashResults.filter(_.hashKey.checkForHitOrMiss)

    if (filtered.nonEmpty) {
      val hashes = CallCacheHashes(filtered)
      val subsets = hashes.hashes.grouped(100)
      subsets foreach { subset =>
        callCacheReadActor ! CacheLookupRequest(CallCacheHashes(subset))
      }
    } else ()
  }

  def updateStateDataWithNewHashResultsAndTransition(hashResults: Set[HashResult]) = {
    if (activity.writeToCache) {
      val newData = stateData.withNewKnownHashes(hashResults)
      if (newData.isDefinitelyCacheHitOrMiss) {
        log.info("New hash results, hit or miss already known (none are cache-checked. Checking if we're done...)")
        checkWhetherAllHashesAreKnownAndTransition(newData)
      } else {
        stay using newData
      }
    } else {
      stay using stateData
    }
  }
}

object EngineJobHashingActor {

  def props(receiver: ActorRef,
            jobDescriptor: BackendJobDescriptor,
            initializationData: Option[BackendInitializationData],
            fileHashingActor: ActorRef,
            callCacheReadActor: ActorRef,
            runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
            backendName: String,
            activity: CallCachingActivity): Props = Props(new EngineJobHashingActor(
      receiver = receiver,
      jobDescriptor = jobDescriptor,
      initializationData = initializationData,
      fileHashingActor = fileHashingActor,
      callCacheReadActor = callCacheReadActor,
      runtimeAttributeDefinitions = runtimeAttributeDefinitions,
      backendName = backendName,
      activity = activity)).withDispatcher(EngineDispatcher)

  private[callcaching] case class EJHAInitialHashingResults(hashes: Set[HashResult]) extends SuccessfulHashResultMessage
  private[callcaching] case object CheckWhetherAllHashesAreKnown

  sealed trait EJHAState
  case object DeterminingHitOrMiss extends EJHAState
  case object GeneratingAllHashes extends EJHAState

  sealed trait EJHAResponse
  case class CacheHit(cacheResultIds: NonEmptyList[CallCachingEntryId]) extends EJHAResponse
  case object CacheMiss extends EJHAResponse
  case class HashError(t: Throwable) extends EJHAResponse {
    override def toString = s"HashError(${t.getMessage})"
  }
  case class CallCacheHashes(hashes: Set[HashResult]) extends EJHAResponse
  object UnspecifiedRuntimeAttributeHashValue extends HashValue("N/A")

  implicit class StringMd5er(unhashedString: String) {
    def md5HashValue: HashValue = {
      val hashBytes = java.security.MessageDigest.getInstance("MD5").digest(unhashedString.getBytes)
      HashValue(javax.xml.bind.DatatypeConverter.printHexBinary(hashBytes))
    }
  }
}

/**
  * Transient data for the EJHA.
  *
  * @param possibleCacheResults The set of cache results which have matched all currently tried hashes
  * @param remainingCacheChecks The set of hash keys which have not yet had their cache results fetched
  * @param hashesKnown The set of all hashes calculated so far (including initial hashes)
  * @param remainingHashesNeeded The set of hashes which are still needed for writing to the database
  */
private[callcaching] case class EJHAData(possibleCacheResults: Option[Set[CallCachingEntryId]],
                                         remainingCacheChecks: Set[HashKey],
                                         hashesKnown: Set[HashResult],
                                         remainingHashesNeeded: Set[HashKey]) {
  // Manipulators
  def intersectCacheResults(newCacheResults: CacheResultMatchesForHashes) = {
    val newIds = newCacheResults.cacheResultIds
    val intersectedIds = possibleCacheResults match {
      case None => newIds
      case Some(currentCacheResults) => currentCacheResults.intersect(newIds)
    }
    this.copy(
      possibleCacheResults = Option(intersectedIds),
      remainingCacheChecks = remainingCacheChecks.diff(newCacheResults.hashResults.map(_.hashKey))
    )
  }
  def withNewKnownHashes(hashResults: Set[HashResult]) = {
    this.copy(
      hashesKnown = hashesKnown ++ hashResults,
      remainingHashesNeeded = remainingHashesNeeded.diff(hashResults.map(_.hashKey)))
  }

  // Queries
  def allHashesKnown = remainingHashesNeeded.isEmpty
  def allCacheResultsIntersected = remainingCacheChecks.isEmpty
  def cacheHit = if (allCacheResultsIntersected) possibleCacheResults flatMap { _.headOption } else None
  def cacheHits = if (allCacheResultsIntersected) possibleCacheResults else None
  def isDefinitelyCacheHit = cacheHits.isDefined
  def isDefinitelyCacheMiss = possibleCacheResults.exists(_.isEmpty)
  def isDefinitelyCacheHitOrMiss = isDefinitelyCacheHit || isDefinitelyCacheMiss
}

private[callcaching] object EJHAData {
  def apply(hashesNeeded: Set[HashKey], activity: CallCachingActivity): EJHAData = EJHAData(
    None,
    if (activity.readFromCache) hashesNeeded.filter(_.checkForHitOrMiss) else Set.empty,
    Set.empty,
    if (activity.writeToCache) hashesNeeded else Set.empty)
}
