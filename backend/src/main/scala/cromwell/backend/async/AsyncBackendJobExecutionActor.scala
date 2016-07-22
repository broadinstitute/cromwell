package cromwell.backend.async

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse, _}
import cromwell.backend.async.AsyncBackendJobExecutionActor._
import cromwell.core.CromwellFatalException
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.services.metadata.MetadataService
import MetadataService.MetadataServiceResponse

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success}

object AsyncBackendJobExecutionActor {

  sealed trait AsyncBackendJobExecutionActorMessage
  private final case class IssuePollRequest(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage
  private final case class PollResponseReceived(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage
  private final case class FailAndStop(reason: Throwable) extends AsyncBackendJobExecutionActorMessage
  private final case class Finish(executionHandle: ExecutionHandle) extends AsyncBackendJobExecutionActorMessage

  trait JobId

  sealed trait ExecutionMode extends AsyncBackendJobExecutionActorMessage {
    def jobId: Option[JobId] = None
  }
  case object Execute extends ExecutionMode
  final case class Recover(recoveryId: JobId) extends ExecutionMode {
    override def jobId = Option(recoveryId)
  }
  final case class UseCachedCall(cachedBackendCall: BackendJobDescriptor) extends ExecutionMode
}


trait AsyncBackendJobExecutionActor { this: Actor with ActorLogging =>

  def retryable: Boolean

  private def withRetry[A](work: () => Future[A], backoff: SimpleExponentialBackoff): Future[A] = {
    def isFatal(t: Throwable) = t.isInstanceOf[CromwellFatalException]

    Retry.withRetry(work, isTransient = !isFatal(_), isFatal = isFatal, backoff = backoff)(context.system)
  }

  private def robustExecuteOrRecover(mode: ExecutionMode) = {
    withRetry(() => executeOrRecover(mode), executeOrRecoverBackoff) onComplete {
      case Success(h) => self ! IssuePollRequest(h)
      case Failure(t) => self ! FailAndStop(t)
    }
  }

  def pollBackoff: SimpleExponentialBackoff

  def executeOrRecoverBackoff: SimpleExponentialBackoff

  private def robustPoll(handle: ExecutionHandle) = {
    withRetry(() => poll(handle), pollBackoff) onComplete {
      case Success(h) => self ! PollResponseReceived(h)
      case Failure(t) => self ! FailAndStop(t)
    }
  }

  private def failAndStop(t: Throwable) = {
    val responseBuilder = if (retryable) FailedRetryableResponse else FailedNonRetryableResponse
    completionPromise.success(responseBuilder.apply(jobDescriptor.key, t, None))
    context.stop(self)
  }

  def receive: Receive = {
    case mode: ExecutionMode => robustExecuteOrRecover(mode)
    case IssuePollRequest(handle) => robustPoll(handle)
    case PollResponseReceived(handle) if handle.isDone => self ! Finish(handle)
    case PollResponseReceived(handle) =>
      context.system.scheduler.scheduleOnce(pollBackoff.backoffMillis.millis, self, IssuePollRequest(handle))
    case Finish(SuccessfulExecutionHandle(outputs, returnCode, hash, resultsClonedFrom)) =>
      completionPromise.success(SucceededResponse(jobDescriptor.key, Some(returnCode), outputs))
      context.stop(self)
    case Finish(FailedNonRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(FailedNonRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(FailedRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(FailedRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(cromwell.backend.async.AbortedExecutionHandle) =>
      completionPromise.success(AbortedResponse(jobDescriptor.key))
      context.stop(self)
    case FailAndStop(t) => failAndStop(t)
    case response: MetadataServiceResponse => handleMetadataServiceResponse(sender(), response)
    case badMessage => log.error(s"Unexpected message $badMessage.")
  }

  /**
    * Handles metadata service responses, with a default implementation that ignores all successes and failures.
    *
    * Any AsyncBackendJobExecutionActor that happens to use the serviceRegistryActor will have ack messages returning.
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
