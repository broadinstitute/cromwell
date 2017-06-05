package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Status}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching.HashResult
import cromwell.database.sql.joins.CallCachingDiffJoin
import cromwell.database.sql.tables.CallCachingEntry
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
      case CallCacheDiffRequest(parameters) => cache.callCacheDiff(parameters.callA, parameters.callB) map {
        case CallCachingDiffJoin(cacheEntryA, cacheEntryB, diff) =>
          val callATag = List(cacheEntryA.workflowExecutionUuid, cacheEntryA.callFullyQualifiedName, cacheEntryA.jobIndex).mkString(":")
          val callBTag = List(cacheEntryB.workflowExecutionUuid, cacheEntryB.callFullyQualifiedName, cacheEntryB.jobIndex).mkString(":")
          
          val formatted = diff map {
            // If both are Some we assume the keys are the same, as they should be
            case (Some((keyA, valueA)), Some((_, valueB))) =>
              Map(
                keyA -> Map(
                  callATag -> Option(valueA),
                  callBTag -> Option(valueB)
                )
              )
            case (Some((keyA, valueA)), None) =>
              Map(
                keyA -> Map(
                  callATag -> Option(valueA),
                  callBTag -> None
                )
              )
            case (None, Some((keyB, valueB))) =>
              Map(
                keyB -> Map(
                  callATag -> None,
                  callBTag -> Option(valueB)
                )
              )
            // Should not be possible...
            case (None, None) =>
              Map.empty[String, Map[String, Option[String]]]
          }
          
          CallCachingDiff(cacheEntryA, cacheEntryB, formatted)
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
  case class CacheLookupRequest(aggregatedCallHashes: AggregatedCallHashes, cacheHitNumber: Int) extends CallCacheReadActorRequest
  case class HasMatchingInitialHashLookup(aggregatedTaskHash: String) extends CallCacheReadActorRequest
  case class HasMatchingInputFilesHashLookup(fileHashes: NonEmptyList[HashResult]) extends CallCacheReadActorRequest
  case class CallCacheDiffRequest(queryParameter: CallCacheDiffQueryParameter) extends CallCacheReadActorRequest
  
  sealed trait CallCacheReadActorResponse
  // Responses on whether or not there is at least one matching entry (can for initial matches of file matches)
  case object HasMatchingEntries extends CallCacheReadActorResponse
  case object NoMatchingEntries extends CallCacheReadActorResponse

  // Responses when asking for the next cache hit
  case class CacheLookupNextHit(hit: CallCachingEntryId) extends CallCacheReadActorResponse
  case object CacheLookupNoHit extends CallCacheReadActorResponse
  
  // Responses for call cache diff
  case class CallCachingDiff(cacheEntryA: CallCachingEntry,
                             cacheEntryB: CallCachingEntry,
                             hashDifferential: Seq[Map[String, Map[String, Option[String]]]]
                            ) extends CallCacheReadActorResponse
  
  // Failure Response
  case class CacheResultLookupFailure(reason: Throwable) extends CallCacheReadActorResponse
}
