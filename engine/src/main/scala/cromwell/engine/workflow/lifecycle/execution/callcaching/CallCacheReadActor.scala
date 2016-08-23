package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.ExecutionContext
import akka.pattern.pipe
import cromwell.core.callcaching.HashResult
import cromwell.database.sql.MetaInfoId
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._

/**
  * Queues up work sent to it because its receive is non-blocking.
  *
  * Would be nice if instead there was a pull- rather than push-based mailbox but I can't find one...
  */
class CallCacheReadActor(cache: CallCache) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  var requestQueue: List[RequestTuple] = List.empty
  var currentRequester: Option[ActorRef] = None

  override def receive: Receive = {
    case CacheLookupRequest(callCacheHashes) =>
      receiveNewRequest(callCacheHashes)
    case r: CallCacheReadActorResponse =>
      currentRequester foreach { _ ! r }
      cycleRequestQueue()
    case other =>
      log.error("Unexpected message type to CallCacheReadActor: " + other.getClass.getSimpleName)
  }

  private def runRequest(callCacheHashes: CallCacheHashes) = {
    val response = cache.fetchMetaInfoIdsMatchingHashes(callCacheHashes) map {
      CacheResultMatchesForHashes(callCacheHashes.hashes, _)
    } recover {
      case t => CacheResultLookupFailure(t)
    }

    response.pipeTo(self)
  }

  private def cycleRequestQueue() = requestQueue match {
    case RequestTuple(replyTo, request) :: tail =>
      currentRequester = Option(replyTo)
      requestQueue = tail
      runRequest(request)
    case Nil =>
      currentRequester = None
  }

  private def receiveNewRequest(callCacheHashes: CallCacheHashes) = currentRequester match {
    case Some(x) => requestQueue :+= RequestTuple(sender, callCacheHashes)
    case None =>
      currentRequester = Option(sender)
      runRequest(callCacheHashes)
  }
}

object CallCacheReadActor {
  def props(callCache: CallCache): Props = Props(new CallCacheReadActor(callCache))

  private[CallCacheReadActor] case class RequestTuple(requester: ActorRef, hashes: CallCacheHashes)

  case class CacheLookupRequest(callCacheHashes: CallCacheHashes)

  sealed trait CallCacheReadActorResponse
  case class CacheResultMatchesForHashes(hashResults: Set[HashResult], cacheResultIds: Set[MetaInfoId]) extends CallCacheReadActorResponse
  case class CacheResultLookupFailure(reason: Throwable) extends CallCacheReadActorResponse
}
