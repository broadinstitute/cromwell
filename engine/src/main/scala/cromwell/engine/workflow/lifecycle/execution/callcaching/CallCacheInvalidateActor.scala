package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.Dispatcher.EngineDispatcher

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

class CallCacheInvalidateActor(callCache: CallCache, cacheId: CallCachingEntryId) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  def receiver = context.parent

  callCache.invalidate(cacheId) onComplete {
    case Success(_) =>
      receiver ! CallCacheInvalidatedSuccess
      context.stop(self)
    case Failure(t) =>
      receiver ! CallCacheInvalidatedFailure(t)
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
case object CallCacheInvalidatedSuccess extends CallCacheInvalidatedResponse
case object CallCacheInvalidationUnnecessary extends CallCacheInvalidatedResponse
case class CallCacheInvalidatedFailure(t: Throwable) extends CallCacheInvalidatedResponse