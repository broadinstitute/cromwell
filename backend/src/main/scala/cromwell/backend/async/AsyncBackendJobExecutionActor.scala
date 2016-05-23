package cromwell.backend.async

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse, _}
import cromwell.backend.async.AsyncBackendJobExecutionActor._
import cromwell.core.CromwellFatalException
import cromwell.core.retry.Backoff
import cromwell.services.MetadataServiceActor.MetadataServiceResponse

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success}
import cromwell.services.MetadataServiceActor.MetadataServiceResponse


object AsyncBackendJobExecutionActor {

  sealed trait AsyncBackendJobExecutionActorMessage
  final case class IssuePollRequest(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage
  final case class PollResponseReceived(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage
  final case class Finish(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage

  sealed trait ExecutionMode extends AsyncBackendJobExecutionActorMessage
  case object Execute extends ExecutionMode
  final case class Resume(executionInfos: Map[String, Option[String]]) extends ExecutionMode
  final case class UseCachedCall(cachedBackendCall: BackendJobDescriptor) extends ExecutionMode
}


trait AsyncBackendJobExecutionActor { this: Actor with ActorLogging =>

  /**
    * Schedule work according to the schedule of the `backoff`.
    */
  protected def scheduleWork(work: => Unit): Unit = {
    val interval = backoff.backoffMillis.millis
    context.system.scheduler.scheduleOnce(interval) {
      work
    }
  }

  def backoff: Backoff

  def retryable: Boolean

  /**
    * If the `work` `Future` completes successfully, perform the `onSuccess` work, otherwise schedule
    * the execution of the `onFailure` work using an exponential backoff.
    */
  def withRetry(work: Future[ExecutionHandle], onSuccess: ExecutionHandle => Unit, onFailure: => Unit): Unit = {
    work onComplete {
      case Success(s) => onSuccess(s)
      case Failure(e: CromwellFatalException) =>
        log.error(e.getMessage, e)
        val responseBuilder = if (retryable) FailedRetryableResponse else FailedNonRetryableResponse
        completionPromise.success(responseBuilder.apply(jobDescriptor.key, e, None))
        context.stop(self)
      case Failure(e: Exception) =>
        log.error(e.getMessage, e)
        scheduleWork(onFailure)
      case Failure(throwable) =>
        // This is a catch-all for a JVM-ending kind of exception, which is why we throw the exception
        log.error(throwable.getMessage, throwable)
        throw throwable
    }
  }

  def receive: Receive = {
    case mode: ExecutionMode =>
      withRetry(executeOrRecover(mode),
        onSuccess = self ! IssuePollRequest(_),
        onFailure = self ! mode
      )

    case IssuePollRequest(handle) =>
      withRetry(poll(handle),
        onSuccess = self ! PollResponseReceived(_),
        onFailure = self ! IssuePollRequest(handle)
      )
    case PollResponseReceived(handle) if handle.isDone => self ! Finish(handle)
    case PollResponseReceived(handle) => scheduleWork(self ! IssuePollRequest(handle))
    case Finish(SuccessfulExecutionHandle(outputs, returnCode, hash, resultsClonedFrom)) =>
      completionPromise.success(SucceededResponse(jobDescriptor.key, Some(returnCode), outputs))
      context.stop(self)
    case Finish(FailedNonRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(FailedNonRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(FailedRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(FailedRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(cromwell.backend.async.AbortedExecutionHandle) => ???

    case response: MetadataServiceResponse => handleMetadataServiceResponse(sender(), response)

    case badMessage => log.error(s"Unexpected message $badMessage.")
  }

  /**
    * Handles metadata service responses, with a default implementation that ignores all successes and failures.
    *
    * Any AsyncBackendJobExecutionActor that happens to mix in ServiceRegistryClient will have ack messages returning.
    * One may optionally handle the ack responses here, or use the default implementation which is to ignore the ack.
    * Sub classes may choose to resend the metadata based on the success or failure response.
    *
    * @param response The response from metadata service, possibly a failure to store the metadata due to a network
    *                 hiccup etc.
    */
  protected def handleMetadataServiceResponse(sentBy: ActorRef, response: MetadataServiceResponse): Unit = {}

  /**
    * Update the ExecutionHandle
    */
  def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def completionPromise: Promise[BackendJobExecutionResponse]

  def jobDescriptor: BackendJobDescriptor

  protected implicit def ec: ExecutionContext
}
