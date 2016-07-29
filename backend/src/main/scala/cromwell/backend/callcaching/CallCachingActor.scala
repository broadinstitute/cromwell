package cromwell.backend.callcaching

import akka.actor.{LoggingFSM, Props}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend.callcaching.CallCachingActor._
import cromwell.core._
import wdl4s.RuntimeAttributes

import scala.concurrent.ExecutionContext

/**
  * The call caching actor is a total maverick!
  *
  * It goes off on its own as soon as it is created and decides whether the JobDescriptor it has received is a cache hit.
  */
case class CallCachingActor
(
  backendJobDescriptor: BackendJobDescriptor,
  hashFileFunction: String => String,
  hashRuntimeAttributesFunction: RuntimeAttributes => String,
  readFromCache: Boolean,
  writeToCache: Boolean
) extends LoggingFSM[CCAState, CCAData] {

  implicit val ec: ExecutionContext = context.dispatcher

  private val anyWorkRequired = readFromCache || writeToCache

  if (anyWorkRequired) calculateHash()

  startWith(
    if (anyWorkRequired) DeterminingCacheHit else NeverDidAnything,
    NoData
  )

  when(DeterminingCacheHit) {
    case Event(CheckCache, _) =>
      sender ! CallCacheCheckInProgress
      stay()
    case Event(CacheHitInternal(cacheResultId), _) if readFromCache =>
      readEntryFromDatabase(cacheResultId)
      goto(FetchingJobResult)
    case Event(CacheHitInternal(_), _) if !readFromCache =>
      goto(UnreadCacheHit)
    case Event(CacheMissInternal(hashes), _) if writeToCache =>
      goto(AwaitingJobResultAndAllHashes) using WithHashes(hashes)
    case Event(CacheMissInternal(_), _) if !writeToCache =>
      goto(UnupdatedCacheMiss)
  }

  when(AwaitingJobResultAndAllHashes) {
    case Event(CheckCache, _) =>
      sender ! CallCacheMiss
      stay()
    case Event(WriteCache(returnCode, jobOutputs), WithHashes(hashes)) =>
      writeEntryToDatabase(returnCode, jobOutputs)
      sender ! CallCacheWriteInProgress
      goto(WritingCacheEntry) using WithJobResult(returnCode, jobOutputs)
  }

  when(WritingCacheEntry) {
    case Event(CheckCache, _) =>
      sender ! CallCacheMiss
      stay()
    case Event(CheckWriteComplete, _) =>
      sender ! CallCacheWriteInProgress
      stay()
    case Event(WriteCache(_, _), _) =>
      // For the sake of idempotency, respond to this message too:
      sender ! CallCacheWriteInProgress
      stay()
    case Event(CallCacheWriteSuccessInternal, oldData @ WithJobResult(_, _)) =>
      goto(NewCacheEntryWritten) using oldData
  }

  when(FetchingJobResult) {
    case Event(CheckCache, _) =>
      sender ! CallCacheCheckInProgress
      stay()
    case Event(CacheHitResultInternal(returnCode, jobOutputs), WithHashes(hashes)) =>
      goto(CacheLookupSuccess) using WithJobResult(returnCode, jobOutputs)
  }

  // Terminal States:

  when(NewCacheEntryWritten) {
    case Event(CheckCache, WithJobResult(returnCode, jobOutputs)) =>
      sender ! CallCacheMiss
      stay()
    case Event(WriteCache, _) =>
      sender ! CallCacheWriteSuccess
      stay()
  }

  when(UnreadCacheHit) {
    case Event(CheckCache, _) =>
      sender ! CallCacheMiss // Let's just pretend we missed.
      stay()
    case Event(WriteCache, _) =>
      sender ! CallCacheWriteSuccess // Pretend it was written. As far as we are concerned, there is an existing call cache entry for this call.
      stay()
  }

  when(UnupdatedCacheMiss) {
    case Event(CheckCache, _) =>
      sender ! CallCacheMiss
      stay()
    case Event(WriteCache, _) =>
      sender ! CallCacheWriteSuccess // We probably shouldn't receive this but if we do we didn't fail... So the client should be happy.
      stay()
  }

  when(CacheLookupSuccess) {
    case Event(CheckCache, WithJobResult(returnCode, jobOutputs)) =>
      sender ! CallCacheHit(SucceededResponse(backendJobDescriptor.key, returnCode, jobOutputs))
      stay()
    case Event(WriteCache, _) =>
      sender ! CallCacheWriteSuccess // We didn't write it, but we didn't fail either. The client should be happy.
      stay()
  }

  // If you never did anything, everything was a no-op.
  when(NeverDidAnything) {
    case Event(CheckCache, _) =>
      sender ! CallCacheMiss
      stay()
    case Event(WriteCache(_, _), _) =>
      sender ! CallCacheWriteSuccess
      stay()
  }

  private def calculateHash(): Unit = {
    // TODO: Do the work and then:
    // self ! CacheHitInternal(cacheResultId)
    // or collate all the hashes and then:
    self ! CacheMissInternal(Map.empty)
  }

  private def writeEntryToDatabase(returnCode: Option[Int], jobOutputs: JobOutputs): Unit = {
    // TODO: Do the work and then:
    self ! CallCacheWriteSuccessInternal
  }

  private def readEntryFromDatabase(cacheResultId: Int): Unit = {
    // TODO: Do the work and then:
    // self ! CacheHitResultInternal(returnCode, jobOutputs)
    ???
  }
}

object CallCachingActor {

  def props
  (
    backendJobDescriptor: BackendJobDescriptor,
    hashFileFunction: String => String,
    hashRuntimeAttributesFunction: RuntimeAttributes => String,
    readFromCache: Boolean,
    writeToCache: Boolean
  ): Props = Props(CallCachingActor(backendJobDescriptor, hashFileFunction, hashRuntimeAttributesFunction, readFromCache = readFromCache, writeToCache = writeToCache))

  sealed trait CCAState
  case object DeterminingCacheHit extends CCAState // We don't know if this is a cache hit yet
  case object AwaitingJobResultAndAllHashes extends CCAState // We know this wasn't a cache hit. Let's wait for some results
  case object WritingCacheEntry extends CCAState // We know this wasn't a cache hit. Let's update the DB with the new results
  case object FetchingJobResult extends CCAState // We know this was a cache hit. Reading the result from the DB

  sealed trait CCATerminalState extends CCAState
  case object NewCacheEntryWritten extends CCATerminalState // We have the hashes, we have the job results. The DB is up to date
  case object UnreadCacheHit extends CCATerminalState // We have the hashes, we had a cache hit, but we aren't reading it
  case object UnupdatedCacheMiss extends CCATerminalState // We have the hashes, we didn't hit the cache, and we won't update it
  case object CacheLookupSuccess extends CCATerminalState // We have the hashes, we hit the cache, and we have retrieved the results
  case object NeverDidAnything extends CCATerminalState // We did nothing, absolutely nothing that day. And I say, what the hell am I doing drinking in LA at 26?

  sealed trait CCAData
  case object NoData extends CCAData
  case class WithHashes(hashes: Map[String, String]) extends CCAData
  case class WithJobResult(returnCode: Option[Int], jobOutputs: JobOutputs) extends CCAData

  sealed trait CCACommand
  case object CheckCache extends CCACommand
  case class WriteCache(returnCode: Option[Int], jobOutputs: JobOutputs) extends CCACommand
  case object CheckWriteComplete extends CCACommand

  sealed trait CCACheckResponse
  case object CallCacheCheckInProgress extends CCACheckResponse
  case class CallCacheHit(backendJobExecutionResponse: BackendJobExecutionResponse) extends CCACheckResponse
  case object CallCacheMiss extends CCACheckResponse

  sealed trait CCAWriteResponse
  case object CallCacheWriteInProgress extends CCAWriteResponse
  case object CallCacheWriteSuccess extends CCAWriteResponse
  case class CallCacheWriteFailure(reason: Throwable) extends CCAWriteResponse

  private[callcaching] sealed trait CCAInternalMessage
  private[callcaching] case class CacheMissInternal(hashes: Map[String, String]) extends CCAInternalMessage
  private[callcaching] case class CacheHitInternal(cacheEntryId: Int) extends CCAInternalMessage
  private[callcaching] case class CacheHitResultInternal(returnCode: Option[Int], jobOutputs: JobOutputs) extends CCAInternalMessage
  private[callcaching] case object CallCacheWriteSuccessInternal extends CCAInternalMessage
  private[callcaching] case class CallCacheWriteFailureInternal(reason: Throwable) extends CCAInternalMessage
}