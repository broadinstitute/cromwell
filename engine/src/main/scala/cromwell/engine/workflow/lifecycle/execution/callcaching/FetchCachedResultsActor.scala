package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.Simpletons._
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object FetchCachedResultsActor {
  def props(cacheHit: CacheHit, replyTo: ActorRef, callCache: CallCache): Props =
    Props(new FetchCachedResultsActor(cacheHit, replyTo, callCache))

  sealed trait CachedResultResponse
  case class CachedOutputLookupFailed(metaInfoId: MetaInfoId, failure: Throwable) extends CachedResultResponse
  case class CachedOutputLookupSucceeded(simpletons: Seq[WdlValueSimpleton], callOutputFiles: Map[String,String],
                                         returnCode: Option[Int], cacheHit: CacheHit, cacheHitDetails: String) extends CachedResultResponse
}


class FetchCachedResultsActor(cacheHit: CacheHit, replyTo: ActorRef, callCache: CallCache)
  extends Actor with ActorLogging {

  {
    implicit val ec: ExecutionContext = context.dispatcher
    val cacheResultId = cacheHit.cacheResultId

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
                                              cacheHit, sourceCacheDetails)
      case Success(None) =>
        val reason = new RuntimeException(s"Cache hit vanished between discovery and retrieval: $cacheResultId")
        replyTo ! CachedOutputLookupFailed(cacheResultId, reason)
      case Failure(t) => replyTo ! CachedOutputLookupFailed(cacheResultId, t)
    }
  }

  override def receive: Receive = Actor.emptyBehavior
}
