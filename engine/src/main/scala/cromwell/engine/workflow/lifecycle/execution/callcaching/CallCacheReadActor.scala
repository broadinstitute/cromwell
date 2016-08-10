package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.database.CromwellDatabase
import cromwell.engine.workflow.lifecycle.execution.{CacheResultLookupFailure, CacheResultMatchesForHashes}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object CallCacheReadActor {
  def props(callCacheHashes: CallCacheHashes): Props = {
    Props(new CallCacheReadActor(callCacheHashes))
  }
}

class CallCacheReadActor(callCacheHashes: CallCacheHashes) extends Actor with ActorLogging {
  {
    implicit val ec: ExecutionContext = context.dispatcher

    val replyTo = context.parent
    val cache = new CallCache(CromwellDatabase.databaseInterface)
    cache.fetchMetaInfoIdsMatchingHashes(callCacheHashes) onComplete {
      case Success(s) => replyTo ! CacheResultMatchesForHashes(callCacheHashes.hashes, s)
      case Failure(t) => replyTo ! CacheResultLookupFailure(t)
    }
  }

  override def receive: Receive = {
    case any => log.error("Unexpected message to CallCacheReadActor: " + any)
  }
}
