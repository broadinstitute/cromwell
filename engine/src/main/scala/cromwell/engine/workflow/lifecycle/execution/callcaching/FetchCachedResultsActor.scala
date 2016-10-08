package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.Simpletons._
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object FetchCachedResultsActor {
  def props(metaInfoId: MetaInfoId, replyTo: ActorRef, callCache: CallCache): Props =
    Props(new FetchCachedResultsActor(metaInfoId, replyTo, callCache))

  sealed trait CachedResultResponse
  case class CachedOutputLookupFailed(metaInfoId: MetaInfoId, failure: Throwable) extends CachedResultResponse
  case class CachedOutputLookupSucceeded(simpletons: Seq[WdlValueSimpleton], callOutputFiles: Map[String,String],
                                         returnCode: Option[Int], cacheHit: MetaInfoId, cacheHitDetails: String) extends CachedResultResponse
}


class FetchCachedResultsActor(cacheResultId: MetaInfoId, replyTo: ActorRef, callCache: CallCache)
  extends Actor with ActorLogging {

  {
    implicit val ec: ExecutionContext = context.dispatcher

    callCache.fetchCachedResult(cacheResultId) onComplete {
      case Success(Some(result)) =>
        val simpletons = result.callCachingSimpletonEntries map toSimpleton
        val jobDetritusFiles = result.callCachingDetritusEntries map { jobDetritusEntry =>
          jobDetritusEntry.detritusKey -> jobDetritusEntry.detritusValue
        }
        val sourceCacheDetails = Seq(result.callCachingEntry.workflowExecutionUuid,
                                    result.callCachingEntry.callFullyQualifiedName,
                                    result.callCachingEntry.jobIndex).mkString(":")

        replyTo ! CachedOutputLookupSucceeded(simpletons, jobDetritusFiles.toMap,
                                              result.callCachingEntry.returnCode,
                                              cacheResultId, sourceCacheDetails)
      case Success(None) =>
        val reason = new RuntimeException(s"Cache hit vanished between discovery and retrieval: $cacheResultId")
        replyTo ! CachedOutputLookupFailed(cacheResultId, reason)
      case Failure(t) => replyTo ! CachedOutputLookupFailed(cacheResultId, t)
    }
  }

  override def receive: Receive = Actor.emptyBehavior
}
