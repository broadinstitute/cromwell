package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.BackendJobExecutionActor
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

case class CallCacheWriteActor(callCache: CallCache, workflowId: WorkflowId, callCacheHashes: CallCacheHashes, succeededResponse: BackendJobExecutionActor.JobSucceededResponse) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  def receiver = context.parent

  callCache.addToCache(workflowId, callCacheHashes, succeededResponse) onComplete {
    case Success(_) =>
      receiver ! CallCacheWriteSuccess
      context.stop(self)
    case Failure(t) =>
      receiver ! CallCacheWriteFailure(t)
      context.stop(self)
  }

  override def receive = {
    case any => log.error("Unexpected message to CallCacheWriteActor: " + any)
  }
}

object CallCacheWriteActor {
  def props(callCache: CallCache, workflowId: WorkflowId, callCacheHashes: CallCacheHashes, succeededResponse: JobSucceededResponse): Props =
    Props(CallCacheWriteActor(callCache, workflowId, callCacheHashes, succeededResponse)).withDispatcher(EngineDispatcher)
}

case object CallCacheWriteSuccess
case class CallCacheWriteFailure(t: Throwable)
