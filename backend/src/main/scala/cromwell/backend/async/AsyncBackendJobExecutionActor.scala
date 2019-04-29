package cromwell.backend.async


import java.util.concurrent.ExecutionException

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.backend.{BackendJobDescriptor, SlowJobWarning}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.async.AsyncBackendJobExecutionActor._
import cromwell.core.CromwellFatalExceptionMarker
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.services.metadata.MetadataService.MetadataServiceResponse

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
  final case class ReconnectToAbort(recoveryId: JobId) extends ExecutionMode {
    override def jobId = Option(recoveryId)
  }
  final case class Reconnect(recoveryId: JobId) extends ExecutionMode {
    override def jobId = Option(recoveryId)
  }
}

trait AsyncBackendJobExecutionActor { this: Actor with ActorLogging with SlowJobWarning =>

  def dockerImageUsed: Option[String]

  // The scala package object (scala/package.scala) contains a neat list of runtime errors that are always going to be fatal.
  // We also consider any Error as fatal, and include the CromwellFatalExceptionMarker so we can mark our own fatal exceptions.
  def isFatal(throwable: Throwable): Boolean = throwable match {
    case _: Error => true
    case _: RuntimeException => true
    case _: InterruptedException => true
    case _: CromwellFatalExceptionMarker => true
    case e: ExecutionException => Option(e.getCause).exists(isFatal)
    case _ => false
  }

  def isTransient(throwable: Throwable): Boolean = false

  private def withRetry[A](work: () => Future[A], backOff: SimpleExponentialBackoff): Future[A] = {
    Retry.withRetry(work, isTransient = isTransient, isFatal = isFatal, backoff = backOff)(context.system)
  }

  private def robustExecuteOrRecover(mode: ExecutionMode) = {
    withRetry(() => executeOrRecover(mode), executeOrRecoverBackOff) onComplete {
      case Success(h) => self ! IssuePollRequest(h)
      case Failure(t) => self ! FailAndStop(t)
    }
  }

  def pollBackOff: SimpleExponentialBackoff

  def executeOrRecoverBackOff: SimpleExponentialBackoff

  private def robustPoll(handle: ExecutionHandle) = {
    withRetry(() => poll(handle), pollBackOff) onComplete {
      case Success(h) => self ! PollResponseReceived(h)
      case Failure(t) => self ! FailAndStop(t)
    }
  }

  private def failAndStop(t: Throwable) = {
    completionPromise.success(JobFailedNonRetryableResponse(jobDescriptor.key, t, None))
    context.stop(self)
  }

  override def receive: Receive = slowJobWarningReceive orElse {
    case mode: ExecutionMode => robustExecuteOrRecover(mode)
    case IssuePollRequest(handle) => robustPoll(handle)
    case PollResponseReceived(handle) if handle.isDone => self ! Finish(handle)
    case PollResponseReceived(handle) =>
      // This should stash the Cancellable someplace so it can be cancelled once polling is complete.
      // -Ywarn-value-discard
      context.system.scheduler.scheduleOnce(pollBackOff.backoffMillis.millis, self, IssuePollRequest(handle))
      ()
    case Finish(SuccessfulExecutionHandle(outputs, returnCode, jobDetritusFiles, executionEvents, _)) =>
      completionPromise.success(JobSucceededResponse(jobDescriptor.key, Some(returnCode), outputs, Option(jobDetritusFiles), executionEvents, dockerImageUsed, resultGenerationMode = RunOnBackend))
      context.stop(self)
    case Finish(FailedNonRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(JobFailedNonRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(FailedRetryableExecutionHandle(throwable, returnCode)) =>
      completionPromise.success(JobFailedRetryableResponse(jobDescriptor.key, throwable, returnCode))
      context.stop(self)
    case Finish(cromwell.backend.async.AbortedExecutionHandle) =>
      completionPromise.success(JobAbortedResponse(jobDescriptor.key))
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
