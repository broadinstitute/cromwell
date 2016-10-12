package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.callcaching.HashResult
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.{CacheLookupRequest, CacheResultLookupFailure, CacheResultMatchesForHashes}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.util.{Failure, Success, Try}

/**
  * Has a set of responses which it will respond with. If it gets a request for anything that it's not expecting to respond to will generate a failure.
  */
class PredictableCallCacheReadActor(responses: Map[String, Set[CallCachingEntryId]]) extends Actor with ActorLogging {

  var responsesRemaining = responses

  override def receive = {
    case CacheLookupRequest(callCacheHashes: CallCacheHashes) =>
      callCacheHashes.hashes.toList match {
        case Nil =>
          throw new Exception("Should never be looking up 0 hash keys!?")
        case head :: tail =>
          val startSet = toTry(head.hashKey.key, responses.get(head.hashKey.key))
          respond(sender, callCacheHashes.hashes, tail.foldLeft(startSet)(resultLookupFolder))
      }
  }

  private def respond(sndr: ActorRef, hashes: Set[HashResult], result: Try[Set[CallCachingEntryId]]) = result match {
    case Success(cacheMatches) => sndr ! CacheResultMatchesForHashes(hashes, cacheMatches)
    case Failure(t) => sndr ! CacheResultLookupFailure(t)
  }

  private def toTry[A](name: String, option: Option[A]): Try[A] = option match {
    case Some(x) => Success(x)
    case None => Failure(new Exception(s"Error looking up response $name!"))
  }

  private def resultLookupFolder(current: Try[Set[CallCachingEntryId]], next: HashResult): Try[Set[CallCachingEntryId]] = current flatMap { c =>
    val lookedUp = toTry(next.hashKey.key, responses.get(next.hashKey.key))
    lookedUp map { l => c.intersect(l) }
  }
}
