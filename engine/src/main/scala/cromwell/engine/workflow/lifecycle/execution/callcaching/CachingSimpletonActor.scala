package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.JobOutputs
import cromwell.database.CromwellDatabase
import cromwell.database.sql.{CachedResult, MetaInfoId}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CachingSimpletonActor.{CachedOutputLookupSucceeded, CachedOutputLookupFailed, CacheResultResponse}
import wdl4s.TaskOutput

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object CachingSimpletonActor {
    def props(metaInfoId: MetaInfoId, taskOutputs: Seq[TaskOutput]): Props =
        Props(new CachingSimpletonActor(metaInfoId, taskOutputs))

  sealed trait CachedResultResponse
  case class CachedOutputLookupFailed(metaInfoId: MetaInfoId, failure: Throwable) extends CachedResultResponse
  case class CachedOutputLookupSucceeded(jobOutputs: JobOutputs) extends CachedResultResponse
}


class CachingSimpletonActor(metaInfoId: MetaInfoId, taskOutputs: Seq[TaskOutput]) extends Actor with ActorLogging with CromwellDatabase with CacheResultResponse {
  {
    implicit val ec: ExecutionContext = context.dispatcher

    val callCache = new CallCache(CromwellDatabase.databaseInterface)
    val replyTo = context.parent

    callCache.fetchCachedResult(metaInfoId) onComplete {
      case Success(Some(s)) => processJobOutputs(s, taskOutputs)
      case Failure(t) => replyTo ! CachedOutputLookupFailed(metaInfoId, t)
    }

    def processJobOutputs(cachedResult: CachedResult, taskOutputs: Seq[TaskOutput]) = {
      callCache.convertToJobOutputs(cachedResult, taskOutputs) match {
        case (jobOutputs: JobOutputs) => replyTo ! CachedOutputLookupSucceeded(jobOutputs)
        case (t: Throwable) => replyTo ! CachedOutputLookupFailed(metaInfoId, t)
      }
    }
  }

  override def receive: Receive = {
    case any => log.error("Unexpected message to CachingSimpletonActor: " + any)
  }

}