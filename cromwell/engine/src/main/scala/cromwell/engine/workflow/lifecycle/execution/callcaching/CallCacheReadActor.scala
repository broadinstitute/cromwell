package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Status}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.core.callcaching.HashResult
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._

import scala.concurrent.ExecutionContext

/**
  * Queues up work sent to it because its receive is non-blocking.
  *
  * Would be nice if instead there was a pull- rather than push-based mailbox but I can't find one...
  */
class CallCacheReadActor(cache: CallCache) extends Actor with ActorLogging {

  private implicit val ec: ExecutionContext = context.dispatcher

  private var requestQueue: List[RequestTuple] = List.empty
  private var currentRequester: Option[ActorRef] = None

  override def receive: Receive = {
    case request: CallCacheReadActorRequest => receiveNewRequest(request)
    case response: CallCacheReadActorResponse =>
      currentRequester foreach { _ ! response }
      cycleRequestQueue()
    case Status.Failure(f) =>
      currentRequester foreach { _ ! CacheResultLookupFailure(new Exception(s"Call Cache query failure: ${f.getMessage}")) }
      cycleRequestQueue()
    case other =>
      log.error("Unexpected message type to CallCacheReadActor: " + other.getClass.getSimpleName)
  }

  private def runRequest(request: CallCacheReadActorRequest): Unit = {
    val response = request match {
      case HasMatchingInitialHashLookup(initialHash) =>
        cache.hasBaseAggregatedHashMatch(initialHash) map {
          case true => HasMatchingEntries
          case false => NoMatchingEntries
        }
      case HasMatchingInputFilesHashLookup(fileHashes) =>
        cache.hasKeyValuePairHashMatch(fileHashes) map {
          case true => HasMatchingEntries
          case false => NoMatchingEntries
        }
      case CacheLookupRequest(aggregatedCallHashes, cacheHitNumber) =>
        cache.callCachingHitForAggregatedHashes(aggregatedCallHashes, cacheHitNumber) map {
          case Some(nextHit) => CacheLookupNextHit(nextHit)
          case None => CacheLookupNoHit
        }
      case call @ CallCacheEntryForCall(workflowId, jobKey) =>
        import cromwell.core.ExecutionIndex._
        cache.cacheEntryExistsForCall(workflowId.toString, jobKey.call.fullyQualifiedName, jobKey.index.fromIndex) map {
          case true => HasCallCacheEntry(call)
          case false => NoCallCacheEntry(call)
        }
    }

    val recovered = response recover {
      case t => CacheResultLookupFailure(t)
    }

    recovered.pipeTo(self)
    ()
  }

  private def cycleRequestQueue() = requestQueue match {
    case RequestTuple(replyTo, request) :: tail =>
      currentRequester = Option(replyTo)
      requestQueue = tail
      runRequest(request)
    case Nil =>
      currentRequester = None
  }

  private def receiveNewRequest(request: CallCacheReadActorRequest): Unit = currentRequester match {
    case Some(_) => requestQueue :+= RequestTuple(sender, request)
    case None =>
      currentRequester = Option(sender)
      runRequest(request)
  }
}

object CallCacheReadActor {
  def props(callCache: CallCache): Props = Props(new CallCacheReadActor(callCache)).withDispatcher(EngineDispatcher)

  private[CallCacheReadActor] case class RequestTuple(requester: ActorRef, request: CallCacheReadActorRequest)

  object AggregatedCallHashes {
    def apply(baseAggregatedHash: String, inputFilesAggregatedHash: String) = {
      new AggregatedCallHashes(baseAggregatedHash, Option(inputFilesAggregatedHash))
    }
  }
  case class AggregatedCallHashes(baseAggregatedHash: String, inputFilesAggregatedHash: Option[String])

  sealed trait CallCacheReadActorRequest
  final case class CacheLookupRequest(aggregatedCallHashes: AggregatedCallHashes, cacheHitNumber: Int) extends CallCacheReadActorRequest
  final case class HasMatchingInitialHashLookup(aggregatedTaskHash: String) extends CallCacheReadActorRequest
  final case class HasMatchingInputFilesHashLookup(fileHashes: NonEmptyList[HashResult]) extends CallCacheReadActorRequest
  final case class CallCacheEntryForCall(workflowId: WorkflowId, jobKey: BackendJobDescriptorKey) extends CallCacheReadActorRequest

  sealed trait CallCacheReadActorResponse
  // Responses on whether or not there is at least one matching entry (can for initial matches of file matches)
  case object HasMatchingEntries extends CallCacheReadActorResponse
  case object NoMatchingEntries extends CallCacheReadActorResponse

  // Responses when asking for the next cache hit
  final case class CacheLookupNextHit(hit: CallCachingEntryId) extends CallCacheReadActorResponse
  case object CacheLookupNoHit extends CallCacheReadActorResponse
  
  final case class HasCallCacheEntry(call: CallCacheEntryForCall) extends CallCacheReadActorResponse
  final case class NoCallCacheEntry(call: CallCacheEntryForCall) extends CallCacheReadActorResponse
  
  // Failure Response
  case class CacheResultLookupFailure(reason: Throwable) extends CallCacheReadActorResponse
}
