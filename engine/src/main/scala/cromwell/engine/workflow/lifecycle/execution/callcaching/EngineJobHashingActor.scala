package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.backend.callcaching.FileHasherWorkerActor.SingleFileHashRequest
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, RuntimeAttributeDefinition}
import cromwell.database.sql.MetaInfoId
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.DockerHashLookupWorkerActor.{DockerHashLookupCommand, DockerHashLookupKey}
import cromwell.engine.workflow.lifecycle.execution.{CacheResultLookupFailure, CacheResultMatchesForHashes}
import wdl4s.values.WdlFile

/**
  * Over time will emit up to two messages back to its parent:
  *  * (if read enabled): Either a CacheHit(id) or CacheMiss message
  *  * (if write enabled): A CallCacheHashes(hashes) message
  */
case class EngineJobHashingActor(jobDescriptor: BackendJobDescriptor,
                                 initializationData: Option[BackendInitializationData],
                                 fileHasherActor: ActorRef,
                                 dockerHashLookupActor: ActorRef,
                                 runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
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
      updateStateDataWithNewHashResultsAndTransition(newHashResults)
    case Event(newCacheResults: CacheResultMatchesForHashes, _) =>
      checkWhetherHitOrMissIsKnownThenTransition(stateData.intersectCacheResults(newCacheResults))
  }

  when(GeneratingAllHashes) {
    case Event(hashResultMessage: SuccessfulHashResultMessage, _) => checkWhetherAllHashesAreKnownAndTransition(stateData.withNewKnownHashes(hashResultMessage.hashes))
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

  onTransition {
    case fromState -> toState =>
      log.debug("Transitioning from {}({}) to {}({})", fromState, stateData, toState, nextStateData)
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
      case HashFileContents => fileInputSimpletons map { case WdlValueSimpleton(name, x: WdlFile) => SingleFileHashRequest(jobDescriptor.key, HashKey(s"input: File $name"), x, initializationData) }
      case HashFilePath => List.empty
    }

    val dockerLookupHashNeeded = (jobDescriptor.runtimeAttributes.get("docker"), mode.dockerHashingType) match {
      case (Some(dockerName), HashDockerNameAndLookupDockerHash) => Set(DockerHashLookupCommand(dockerName.valueString))
      case _ => Set.empty
    }

    val hashesNeeded: Set[HashKey] = initialHashes.map(_.hashKey) ++ fileContentHashesNeeded.map(_.hashKey) ++ dockerLookupHashNeeded.map(_ => DockerHashLookupKey)

    val initialState = if (mode.readFromCache) DeterminingHitOrMiss else GeneratingAllHashes
    val initialData = EJHAData(jobDescriptor, hashesNeeded, mode)

    startWith(initialState, initialData)

    // Submit the set of initial hashes for checking against the DB:
    self ! EJHAInitialHashingResults(initialHashes)
    // Find the hashes for all input files:
    fileContentHashesNeeded.foreach(fileHasherActor ! _)
    dockerLookupHashNeeded.foreach(dockerHashLookupActor ! _)
  }

  private def calculateInitialHashes(nonFileInputs: Iterable[WdlValueSimpleton], fileInputs: Iterable[WdlValueSimpleton]): Set[HashResult] = {

    val commandTemplateHash = HashResult(HashKey("command template"), jobDescriptor.call.task.commandTemplate.toString.md5HashValue)
    val backendNameHash = HashResult(HashKey("backend name"), backendName.md5HashValue)
    val inputCountHash = HashResult(HashKey("input count"), (nonFileInputs.size + fileInputs.size).toString.md5HashValue)

    val runtimeAttributeHashes = runtimeAttributeDefinitions map { definition => jobDescriptor.runtimeAttributes.get(definition.name) match {
      case Some(wdlValue) => HashResult(HashKey("runtime attribute: " + definition.name, definition.usedInCallCaching), wdlValue.valueString.md5HashValue)
      case None => HashResult(HashKey("runtime attribute: " + definition.name, definition.usedInCallCaching), UnspecifiedRuntimeAttributeHashValue)
    }}

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
    Set(commandTemplateHash, backendNameHash, inputCountHash) ++ runtimeAttributeHashes ++ inputHashResults ++ outputExpressionHashResults ++ inputFilePathHashes
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
  private def findCacheResults(hashResults: Set[HashResult]) = {
    val filtered = hashResults.filter(_.hashKey.checkForHitOrMiss)

    if (filtered.nonEmpty) {
      val hashes = CallCacheHashes(filtered)
      val actorName = s"CallCacheReadActor-${jobDescriptor.workflowDescriptor.id}-${jobDescriptor.key.tag}-lookup_$lookupIndex"
      lookupIndex += 1
      context.actorOf(CallCacheReadActor.props(hashes), actorName)
    } else ()
  }

  def updateStateDataWithNewHashResultsAndTransition(hashResults: Set[HashResult]) = {
    if (mode.writeToCache) {
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

  def props(jobDescriptor: BackendJobDescriptor,
            initializationData: Option[BackendInitializationData],
            fileHasherActor: ActorRef,
            dockerHashLookupActor: ActorRef,
            runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
            backendName: String,
            activity: CallCachingActivity): Props =
    Props(new EngineJobHashingActor(jobDescriptor, initializationData, fileHasherActor, dockerHashLookupActor, runtimeAttributeDefinitions, backendName, activity))

  private[callcaching] case class EJHAInitialHashingResults(hashes: Set[HashResult]) extends SuccessfulHashResultMessage
  private[callcaching] case object CheckWhetherAllHashesAreKnown

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
private[callcaching] case class EJHAData(possibleCacheResults: Option[Set[MetaInfoId]],
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
  def isDefinitelyCacheHit = cacheHit.isDefined
  def isDefinitelyCacheMiss = possibleCacheResults.exists(_.isEmpty)
  def isDefinitelyCacheHitOrMiss = isDefinitelyCacheHit || isDefinitelyCacheMiss
}

object EJHAData {
  def apply(jobDescriptor: BackendJobDescriptor, hashesNeeded: Set[HashKey], activity: CallCachingActivity): EJHAData = EJHAData(
    None,
    if (activity.readFromCache) hashesNeeded.filter(_.checkForHitOrMiss) else Set.empty,
    Set.empty,
    if (activity.writeToCache) hashesNeeded else Set.empty)
}
