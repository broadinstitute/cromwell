package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.database.sql.tables.CallCachingEntry
import cromwell.services.CallCaching.CallCachingEntryId

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

class CallCacheInvalidateActor(callCache: CallCache, cacheId: CallCachingEntryId) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  def receiver = context.parent

  callCache.invalidate(cacheId) onComplete {
    case Success(maybeEntry) =>
      receiver ! CallCacheInvalidatedSuccess(cacheId, maybeEntry)
      context.stop(self)
    case Failure(t) =>
      receiver ! CallCacheInvalidatedFailure(cacheId, t)
      context.stop(self)
  }

  override def receive: Receive = {
    case any => log.error("Unexpected message to InvalidateCallCacheActor: " + any)
  }
}

object CallCacheInvalidateActor {
  def props(callCache: CallCache, cacheId: CallCachingEntryId) = {
    Props(new CallCacheInvalidateActor(callCache: CallCache, cacheId: CallCachingEntryId)).withDispatcher(EngineDispatcher)
  }
}

sealed trait CallCacheInvalidatedResponse
case class CallCacheInvalidatedSuccess(cacheId: CallCachingEntryId, maybeEntry: Option[CallCachingEntry]) extends CallCacheInvalidatedResponse
case object CallCacheInvalidationUnnecessary extends CallCacheInvalidatedResponse
case class CallCacheInvalidatedFailure(cacheId: CallCachingEntryId, t: Throwable) extends CallCacheInvalidatedResponse
