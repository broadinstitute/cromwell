package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobCachingActor._
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.JobOutputs

import scala.concurrent.Future

//1. Unclear if Abort needs to be a part of Caching actor? What should happen when an abort command hits during caching and not execution?
//2. Send back the same kind of responses that the BackendJobExecutionActor would send.


object BackendJobCachingActor {

  sealed trait BackendJobCachingActorCommand extends BackendWorkflowLifecycleActorCommand
  case class CacheJobCommand(cachedJobOutputs: JobOutputs) extends BackendJobCachingActorCommand

  sealed trait BackendJobCachingActorResponse extends BackendWorkflowLifecycleActorCommand
  sealed trait BackendJobExecutionResponse extends BackendJobCachingActorResponse { def jobKey: BackendJobDescriptorKey }
}

trait BackendJobCachingActor extends Actor with ActorLogging with BackendJobLifecycleActor {

  def receive: Receive = LoggingReceive {
    case CacheJobCommand => performActionThenRespond(copyCachedOutputs, onFailure = cachingFailed)
//    case AbortJobCommand =>
//      abort()
//      context.parent ! AbortedResponse(jobDescriptor.key)
//      context.stop(self)
  }

  /**
    * Abort a job that's caching.
    */
  def abort(): Unit = {
    log.warning("Aborts not supported for a job being cached.",
      jobTag, jobDescriptor.key.call.fullyQualifiedName)
  }

  def copyCachedOutputs: Future[BackendJobExecutionResponse]

  private def cachingFailed = (t: Throwable) =>
    FailedNonRetryableResponse(jobKey = jobDescriptor.key, throwable = t, returnCode = None)
}