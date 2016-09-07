package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.Simpletons._
import cromwell.core.JobOutputs
import cromwell.core.simpleton.WdlValueBuilder
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.services.SingletonServicesStore
import wdl4s.TaskOutput

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

object FetchCachedResultsActor {
  def props(cacheHit: CacheHit, taskOutputs: Seq[TaskOutput], replyTo: ActorRef, callCache: CallCache): Props =
    Props(new FetchCachedResultsActor(cacheHit, taskOutputs, replyTo, callCache))

  sealed trait CachedOutputLookupResponse
  case class CachedOutputLookupFailed(metaInfoId: MetaInfoId, failure: Throwable) extends CachedOutputLookupResponse
  case class CachedOutputLookupSucceeded(jobOutputs: JobOutputs, cacheHit: CacheHit) extends CachedOutputLookupResponse
}


class FetchCachedResultsActor(cacheHit: CacheHit, taskOutputs: Seq[TaskOutput], replyTo: ActorRef, callCache: CallCache)
  extends Actor with ActorLogging {
  {
    implicit val ec: ExecutionContext = context.dispatcher
    val cacheResultId = cacheHit.cacheResultId

    callCache.fetchCachedResult(cacheResultId) onComplete {
      case Success(Some(result)) =>
        val simpletons = result.resultSimpletons map toSimpleton
        val jobOutputs = WdlValueBuilder.toJobOutputs(taskOutputs, simpletons)
        replyTo ! CachedOutputLookupSucceeded(jobOutputs, cacheHit)
      case Success(None) =>
        val reason = new RuntimeException(s"Cache hit vanished between discovery and retrieval: $cacheResultId")
        replyTo ! CachedOutputLookupFailed(cacheResultId, reason)
      case Failure(t) => replyTo ! CachedOutputLookupFailed(cacheResultId, t)
    }
  }

  override def receive: Receive = Actor.emptyBehavior
}
