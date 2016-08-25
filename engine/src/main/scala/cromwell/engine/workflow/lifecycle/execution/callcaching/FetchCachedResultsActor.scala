package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.JobOutputs
import cromwell.database.CromwellDatabase
import cromwell.database.sql.{CachedResult, MetaInfoId}
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedResultResponse, CachedOutputLookupSucceeded, CachedOutputLookupFailed, CacheResultResponse}
import wdl4s.TaskOutput

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Try}

object FetchCachedResultsActor {
    def props(metaInfoId: MetaInfoId, taskOutputs: Seq[TaskOutput]): Props =
        Props(new FetchCachedResultsActor(metaInfoId, taskOutputs))

  sealed trait CachedResultResponse
  case class CachedOutputLookupFailed(metaInfoId: MetaInfoId, failure: Throwable) extends CachedResultResponse
  case class CachedOutputLookupSucceeded(jobOutputs: JobOutputs) extends CachedResultResponse
}


class FetchCachedResultsActor(metaInfoId: MetaInfoId, taskOutputs: Seq[TaskOutput]) extends Actor with ActorLogging with CromwellDatabase {
  {
    implicit val ec: ExecutionContext = context.dispatcher

    val callCache = new CallCache(CromwellDatabase.databaseInterface)
    val replyTo = context.parent

    val response: Future[CachedResultResponse] = fetchCachedJobResult map { CachedOutputLookupSucceeded(_) } recover { case t => CachedOutputLookupFailed(metaInfoId, t) }

    def fetchCachedJobResult = callCache.fetchCachedResult(metaInfoId) flatMap {
      case Some(r) => Future.fromTry(Try(callCache.convertToJobOutputs(r, taskOutputs)))
      case None => Future.failed(new IllegalArgumentException(s"Missing Cached Results for MetaInfoId: ${metaInfoId}"))
    }

    replyTo ! response
  }

  override def receive: Receive = {
    case any => log.error(s"Unexpected message to ${this.getClass.getName} " + any)
  }

}