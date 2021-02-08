package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.Simpletons._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.database.sql.SqlConverters._
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.services.CallCaching.CallCachingEntryId

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object FetchCachedResultsActor {
  def props(callCachingEntryId: CallCachingEntryId, replyTo: ActorRef, callCache: CallCache): Props =
    Props(new FetchCachedResultsActor(callCachingEntryId, replyTo, callCache)).withDispatcher(EngineDispatcher)

  sealed trait CachedResultResponse
  case class CachedOutputLookupFailed(callCachingEntryId: CallCachingEntryId, failure: Throwable) extends CachedResultResponse
  case class CachedOutputLookupSucceeded(simpletons: Seq[WomValueSimpleton], callOutputFiles: Map[String,String],
                                         returnCode: Option[Int], cacheHit: CallCachingEntryId, cacheHitDetails: String) extends CachedResultResponse
}


class FetchCachedResultsActor(cacheResultId: CallCachingEntryId, replyTo: ActorRef, callCache: CallCache)
  extends Actor with ActorLogging {

  {
    implicit val ec: ExecutionContext = context.dispatcher

    callCache.fetchCachedResult(cacheResultId) map {
      case Some(result) =>
        val simpletons = result.callCachingSimpletonEntries map toSimpleton
        val jobDetritusFiles = result.callCachingDetritusEntries map { jobDetritusEntry =>
          jobDetritusEntry.detritusKey -> jobDetritusEntry.detritusValue.toRawString
        }

        val sourceCacheDetails = Seq(result.callCachingEntry.workflowExecutionUuid,
          result.callCachingEntry.callFullyQualifiedName,
          result.callCachingEntry.jobIndex).mkString(":")

        CachedOutputLookupSucceeded(simpletons, jobDetritusFiles.toMap,
          result.callCachingEntry.returnCode,
          cacheResultId, sourceCacheDetails)
      case None =>
        val reason = new RuntimeException(s"Cache hit vanished between discovery and retrieval: $cacheResultId")
        CachedOutputLookupFailed(cacheResultId, reason)
    } onComplete {
      case Success(sendMe) => replyTo ! sendMe
      case Failure(t) => replyTo ! CachedOutputLookupFailed(cacheResultId, t)
    }
  }

  override def receive: Receive = Actor.emptyBehavior
}
