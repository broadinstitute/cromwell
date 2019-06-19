package cromwell.engine.workflow.lifecycle.execution.callcaching

import java.security.MessageDigest

import akka.actor.{ActorRef, LoggingFSM, Props, Terminated}
import cats.data.NonEmptyList
import cromwell.backend.standard.callcaching.StandardFileHashingActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, RuntimeAttributeDefinition}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching._
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.CallCacheHashingJobActorData._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheMiss
import javax.xml.bind.DatatypeConverter
import wom.RuntimeAttributesKeys
import wom.types._
import wom.values._


/**
  * Actor responsible for calculating individual as well as aggregated hashes for a job.
  * First calculate the initial hashes (individual and aggregated), and send them to its parent
  * as well as the callCacheReadingJobActor if one was provided.
  * From there, either wait for a NextBatchOfFileHashesRequest message from the callCacheReadingJobActor to trigger the next
  * batch of file hash requests, or send it itself if there is no callCacheReadingJobActor.
  * When we get all the file hashes for one batch, send those hashes to the CCRead actor and wait for the next NextBatchOfFileHashesRequest.
  * If it was the last batch and all the files have been hashed, send all the hashes along with the aggregated file hash.
  * If at any point the callCacheReadingJobActor dies (either because it decided this was a CacheMiss or it found a successful
  * CacheHit), either keep hashing the files if writeToCache is true, or die if it's not.
  * Both this actor and the CCRJA are children of the EJHA. 
  * The CCRJA reports its findings (cache hit / miss) directly to the EJHA that forwards them to the EJEA.
  * This actor does not need nor care about cache hit/miss.
  */
class CallCacheHashingJobActor(jobDescriptor: BackendJobDescriptor,
                               callCacheReadingJobActor: Option[ActorRef],
                               initializationData: Option[BackendInitializationData],
                               runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                               backendNameForCallCachingPurposes: String,
                               fileHashingActorProps: Props,
                               callCachingEligible: CallCachingEligible,
                               callCachingActivity: CallCachingActivity,
                               callCachePathPrefixes: Option[CallCachePathPrefixes]
                              ) extends LoggingFSM[CallCacheHashingJobActorState, CallCacheHashingJobActorData] {

  val fileHashingActor = makeFileHashingActor()

  // Watch the read actor, as it will die when it's done (cache miss or successful cache hit)
  // When that happens we want to either stop if writeToCache is false, or keep going
  callCacheReadingJobActor foreach context.watch

  initializeCCHJA()

  override def preStart(): Unit = {
    if (callCacheReadingJobActor.isEmpty && !callCachingActivity.writeToCache) {
      log.error("Programmer error ! There is no reason to have a hashing actor if both read and write to cache are off")
      context.parent ! CacheMiss
      context stop self
    }
    super.preStart()
  }

  when(WaitingForHashFileRequest) {
    case Event(NextBatchOfFileHashesRequest, data) =>
      data.fileHashRequestsRemaining.headOption match {
        case Some(batch) if batch.nonEmpty =>
          batch foreach { fileHashingActor ! _ }
          goto(HashingFiles)
        case _ =>
          sendToCallCacheReadingJobActor(NoFileHashesResult, data)
          stopAndStay(Option(NoFileHashesResult))
      }
    case Event(Terminated(_), data) if callCachingActivity.writeToCache =>
      self ! NextBatchOfFileHashesRequest
      stay() using data.copy(callCacheReadingJobActor = None)
  }

  when(HashingFiles) {
    case Event(FileHashResponse(result), data) =>
      addFileHash(result, data) match {
        case (newData, Some(_: PartialFileHashingResult)) =>
          self ! NextBatchOfFileHashesRequest
          goto(WaitingForHashFileRequest) using newData
        case (newData, Some(finalResult: FinalFileHashingResult)) =>
          sendToCallCacheReadingJobActor(finalResult, newData)
          stopAndStay(Option(finalResult))
        case (newData, None) =>
          stay() using newData
      }
    case Event(Terminated(_), data) if callCachingActivity.writeToCache =>
      stay() using data.copy(callCacheReadingJobActor = None)
  }

  whenUnhandled {
    case Event(Terminated(_), _) =>
      stopAndStay(None)
    case Event(error: HashingFailedMessage, data) =>
      log.error("""Failed to hash "{}": {}""", error.file, error.reason.getMessage)
      sendToCallCacheReadingJobActor(error, data)
      context.parent ! error
      stopAndStay(None)
  }

  // In its own function so it can be overridden in the test
  private [callcaching] def addFileHash(hashResult: HashResult, data: CallCacheHashingJobActorData) = {
    data.withFileHash(hashResult)
  }

  private def stopAndStay(fileHashResult: Option[FinalFileHashingResult]) = {
    fileHashResult foreach { context.parent ! _ }
    context stop fileHashingActor
    context stop self
    stay()
  }

  private def sendToCallCacheReadingJobActor(message: Any, data: CallCacheHashingJobActorData): Unit = {
    data.callCacheReadingJobActor foreach { _ ! message }
  }

  private def initializeCCHJA(): Unit = {
    import cromwell.core.simpleton.WomValueSimpleton._

    val unqualifiedInputs = jobDescriptor.evaluatedTaskInputs map { case (declaration, value) => declaration.name -> value }

    val inputSimpletons = unqualifiedInputs.simplifyForCaching
    val (fileInputSimpletons, nonFileInputSimpletons) = inputSimpletons partition {
      case WomValueSimpleton(_, _: WomFile) => true
      case _ => false
    }

    val initialHashes = calculateInitialHashes(nonFileInputSimpletons, fileInputSimpletons)

    val fileHashRequests = fileInputSimpletons collect {
      case WomValueSimpleton(name, x: WomFile) => SingleFileHashRequest(jobDescriptor.key, HashKey(true, "input", s"File $name"), x, initializationData)
    }

    val hashingJobActorData = CallCacheHashingJobActorData(fileHashRequests.toList, callCacheReadingJobActor)
    startWith(WaitingForHashFileRequest, hashingJobActorData)

    val aggregatedBaseHash = calculateHashAggregation(initialHashes, MessageDigest.getInstance("MD5"))

    val initialHashingResult = InitialHashingResult(initialHashes, aggregatedBaseHash, callCachePathPrefixes.toList)

    sendToCallCacheReadingJobActor(initialHashingResult, hashingJobActorData)
    context.parent ! initialHashingResult

    // If there is no CCRead actor, we need to send ourselves the next NextBatchOfFileHashesRequest
    if (hashingJobActorData.callCacheReadingJobActor.isEmpty) self ! NextBatchOfFileHashesRequest
  }

  private def calculateInitialHashes(nonFileInputs: Iterable[WomValueSimpleton], fileInputs: Iterable[WomValueSimpleton]): Set[HashResult] = {

    val commandTemplateHash = HashResult(HashKey("command template"), jobDescriptor.taskCall.callable.commandTemplateString(jobDescriptor.evaluatedTaskInputs).md5HashValue)
    val backendNameHash = HashResult(HashKey("backend name"), backendNameForCallCachingPurposes.md5HashValue)
    val inputCountHash = HashResult(HashKey("input count"), (nonFileInputs.size + fileInputs.size).toString.md5HashValue)
    val outputCountHash = HashResult(HashKey("output count"), jobDescriptor.taskCall.callable.outputs.size.toString.md5HashValue)

    val runtimeAttributeHashes = runtimeAttributeDefinitions map { definition => jobDescriptor.runtimeAttributes.get(definition.name) match {
      case Some(_) if definition.name == RuntimeAttributesKeys.DockerKey && callCachingEligible.dockerHash.isDefined =>
        HashResult(HashKey(definition.usedInCallCaching, "runtime attribute", definition.name), callCachingEligible.dockerHash.get.md5HashValue)
      case Some(womValue) => HashResult(HashKey(definition.usedInCallCaching, "runtime attribute", definition.name), womValue.valueString.md5HashValue)
      case None => HashResult(HashKey(definition.usedInCallCaching, "runtime attribute", definition.name), UnspecifiedRuntimeAttributeHashValue)
    }}

    val inputHashResults = nonFileInputs map {
      case WomValueSimpleton(name, value) =>
        val womTypeHashKeyString = value.womType.toHashKeyString
        log.debug("Hashing input expression as {} {}", womTypeHashKeyString, name)
        HashResult(HashKey("input", s"$womTypeHashKeyString $name"),  value.toWomString.md5HashValue)
    }

    val outputExpressionHashResults = jobDescriptor.taskCall.callable.outputs map { output =>
      val womTypeHashKeyString = output.womType.toHashKeyString
      val outputExpressionCacheString = output.expression.cacheString
      log.debug("Hashing output expression type as '{}' and value as '{}'", womTypeHashKeyString, outputExpressionCacheString)
      HashResult(HashKey("output expression", s"$womTypeHashKeyString ${output.name}"), outputExpressionCacheString.md5HashValue)
    }

    // Build these all together for the final set of initial hashes:
    Set(commandTemplateHash, backendNameHash, inputCountHash, outputCountHash) ++ runtimeAttributeHashes ++ inputHashResults ++ outputExpressionHashResults
  }

  private [callcaching] def makeFileHashingActor() = {
    val fileHashingActorName = s"FileHashingActor_for_${jobDescriptor.key.tag}"
    context.actorOf(fileHashingActorProps, fileHashingActorName)
  }
}

object CallCacheHashingJobActor {

  def props(jobDescriptor: BackendJobDescriptor,
            callCacheReadingJobActor: Option[ActorRef],
            initializationData: Option[BackendInitializationData],
            runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
            backendNameForCallCachingPurposes: String,
            fileHashingActorProps: Props,
            callCachingEligible: CallCachingEligible,
            callCachingActivity: CallCachingActivity,
            callCachePathPrefixes: Option[CallCachePathPrefixes]
           ) = Props(new CallCacheHashingJobActor(
    jobDescriptor,
    callCacheReadingJobActor,
    initializationData,
    runtimeAttributeDefinitions,
    backendNameForCallCachingPurposes,
    fileHashingActorProps,
    callCachingEligible,
    callCachingActivity,
    callCachePathPrefixes
  )).withDispatcher(EngineDispatcher)

  sealed trait CallCacheHashingJobActorState
  case object WaitingForHashFileRequest extends CallCacheHashingJobActorState
  case object HashingFiles extends CallCacheHashingJobActorState

  /**
    * Aggregates hash results together in a predictable and reproducible manner.
    * Purposefully takes an Iterable, because it will be sorted appropriately in this method
    * to ensure deterministic result, so the type of collection doesn't matter.
    * Only aggregates hash values for which checkForHitOrMiss is true.
    * If several hash keys are identical, the result of this method is undefined.
    */
  private def calculateHashAggregation(hashes: Iterable[HashResult], messageDigest: MessageDigest) = {
    val sortedHashes = hashes.toList
      .filter(_.hashKey.checkForHitOrMiss)
      .sortBy(_.hashKey.key)
      .map({ case HashResult(hashKey, HashValue(hashValue)) => hashKey.key + hashValue })
      .map(_.getBytes)
    sortedHashes foreach messageDigest.update
    DatatypeConverter.printHexBinary(messageDigest.digest())
  }

  object CallCacheHashingJobActorData {
    // Slick will eventually build a prepared statement with that many parameters. Don't set this too high or it will stackoverflow.
    val BatchSize = 100

    def apply(fileHashRequestsRemaining: List[SingleFileHashRequest], callCacheReadingJobActor: Option[ActorRef]) = {
      new CallCacheHashingJobActorData(fileHashRequestsRemaining.grouped(BatchSize).toList, List.empty, callCacheReadingJobActor)
    }
  }

  final case class CallCacheHashingJobActorData(
                                           fileHashRequestsRemaining: List[List[SingleFileHashRequest]],
                                           fileHashResults: List[HashResult],
                                           callCacheReadingJobActor: Option[ActorRef]
                                         ) {
    private val md5Digest = MessageDigest.getInstance("MD5")

    /**
      * Returns the updated state data along with an optional message to be sent back to CCRJA and parent.
      */
    def withFileHash(hashResult: HashResult): (CallCacheHashingJobActorData, Option[CCHJAFileHashResponse]) = {
      // Add the hash result to the list of known hash results
      val newFileHashResults = hashResult +: fileHashResults

      // Returns a pair of the updated fileHashRequestsRemaining and optionally either a PartialFileHashingResult or a CompleteFileHashingResult
      val (updatedRequestsList, responseMessage) = fileHashRequestsRemaining match {
        case Nil => (List.empty, Option(NoFileHashesResult))
        case lastBatch :: Nil =>
          val updatedBatch = lastBatch.filterNot(_.hashKey == hashResult.hashKey)
          // If we're processing the last batch, and it's now empty, then we're done
          // In that case compute the aggregated hash and send that
          if (updatedBatch.isEmpty) (List.empty, Option(CompleteFileHashingResult(newFileHashResults.toSet, calculateHashAggregation(newFileHashResults, md5Digest))))
          // Otherwise just return the updated batch and no message
          else (List(updatedBatch), None)
        case currentBatch :: otherBatches =>
          val updatedBatch = currentBatch.filterNot(_.hashKey == hashResult.hashKey)
          // If the current batch is empty, we got a partial result, take the first BatchSize of the list
          if (updatedBatch.isEmpty) {
            // hashResult + fileHashResults.take(BatchSize - 1) -> BatchSize elements
            val partialHashes = NonEmptyList.of[HashResult](hashResult, fileHashResults.take(BatchSize - 1): _*)
            (otherBatches, Option(PartialFileHashingResult(partialHashes)))
          }
          // Otherwise just return the updated request list and no message
          else (updatedBatch :: otherBatches, None)
      }

      (this.copy(fileHashRequestsRemaining = updatedRequestsList, fileHashResults = newFileHashResults), responseMessage)
    }
  }

  object UnspecifiedRuntimeAttributeHashValue extends HashValue("N/A")

  sealed trait CCHJARequest
  case object NextBatchOfFileHashesRequest extends CCHJARequest

  sealed trait CCHJAResponse
  case class InitialHashingResult(initialHashes: Set[HashResult], aggregatedBaseHash: String, cacheHitHints: List[CacheHitHint] = List.empty) extends CCHJAResponse

  // File Hashing responses
  sealed trait CCHJAFileHashResponse extends CCHJAResponse
  case class PartialFileHashingResult(initialHashes: NonEmptyList[HashResult]) extends CCHJAFileHashResponse

  sealed trait FinalFileHashingResult extends CCHJAFileHashResponse
  case class CompleteFileHashingResult(fileHashes: Set[HashResult], aggregatedFileHash: String) extends FinalFileHashingResult
  case object NoFileHashesResult extends FinalFileHashingResult

  implicit class StringMd5er(val unhashedString: String) extends AnyVal {
    def md5HashValue: HashValue = {
      val hashBytes = java.security.MessageDigest.getInstance("MD5").digest(unhashedString.getBytes)
      HashValue(javax.xml.bind.DatatypeConverter.printHexBinary(hashBytes))
    }
  }

  implicit class WomTypeHashString(val womType: WomType) extends AnyVal {
    def toHashKeyString: String = {
      womType match {
        case c: WomCompositeType =>
          val fieldTypes = c.typeMap map {
            case (key, value) => s"$key -> ${value.stableName}"
          }
          "CompositeType_digest_" + fieldTypes.mkString("\n").md5Sum
        case a: WomArrayType =>
          s"Array(${a.memberType.toHashKeyString})"
        case o: WomOptionalType =>
          s"Optional(${o.memberType.toHashKeyString})"
        case p: WomPairType =>
          s"Pair(${p.leftType.toHashKeyString},${p.rightType.toHashKeyString})"
        case c: WomCoproductType =>
          val hashStrings = c.types.toList.map(_.toHashKeyString).mkString(",")
          s"Coproduct($hashStrings)"
        case o => o.stableName
      }
    }
  }
}
