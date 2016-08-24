package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse, FailedNonRetryableResponse}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.simpleton.WdlValueSimpleton

import scala.concurrent.Future

object BackendCacheHitCopyingActor {
  final case class CopyOutputsCommand(wdlValueSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String], returnCode: Option[Int])
}

trait BackendCacheHitCopyingActor extends Actor with ActorLogging with BackendJobLifecycleActor {

  def copyCachedOutputs(wdlValueSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String,String], returnCode: Option[Int]): Future[BackendJobExecutionResponse]

  def receive: Receive = LoggingReceive {
    case CopyOutputsCommand(simpletons, jobDetritus, returnCode) =>
      performActionThenRespond(copyCachedOutputs(simpletons, jobDetritus, returnCode), onFailure = cachingFailed, andThen = context stop self)
    case AbortJobCommand =>
      abort()
      context.parent ! AbortedResponse(jobDescriptor.key)
      context stop self
  }

  def abort(): Unit = log.warning("{}: Abort not supported during cache hit copying", jobTag)

  private def cachingFailed(t: Throwable) = {
    FailedNonRetryableResponse(jobKey = jobDescriptor.key, throwable = t, returnCode = None)
  }
}
