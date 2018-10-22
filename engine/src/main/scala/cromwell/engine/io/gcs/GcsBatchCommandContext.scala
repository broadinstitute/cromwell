package cromwell.engine.io.gcs

import java.io.FileNotFoundException

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.client.util.ExponentialBackOff
import com.google.cloud.storage.StorageException
import common.util.Backoff
import cromwell.core.io.SingleFileIoCommand
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.engine.io.IoActor.IoResult
import cromwell.engine.io.gcs.GcsBatchCommandContext.BatchResponse
import cromwell.engine.io.{IoActor, IoCommandContext}
import cromwell.filesystems.gcs.batch.GcsBatchIoCommand
import cromwell.filesystems.gcs.RequesterPaysErrors._

import scala.concurrent.Promise
import scala.concurrent.duration._
import scala.language.postfixOps

object GcsBatchCommandContext {
  def defaultBackoff = SimpleExponentialBackoff(
    new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(1.second.toMillis.toInt)
      .setMultiplier(4)
      .setMaxIntervalMillis(30.seconds.toMillis.toInt)
      .setRandomizationFactor(0.2D)
      .setMaxElapsedTimeMillis(30.minutes.toMillis.toInt)
      .build()
  )
  type BatchResponse = Either[IoResult, GcsBatchCommandContext[_, _]]
}

final case class GcsBatchCommandContext[T, U](request: GcsBatchIoCommand[T, U],
                                              replyTo: ActorRef,
                                              override val clientContext: Option[Any] = None,
                                              backoff: Backoff = GcsBatchCommandContext.defaultBackoff,
                                              currentAttempt: Int = 1,
                                              promise: Promise[BatchResponse] = Promise[BatchResponse]
                                             ) extends IoCommandContext[T] {

  /**
    * None if no retry should be attempted, Some(timeToWaitBeforeNextAttempt) otherwise
    */
  lazy val retryIn = if (currentAttempt >= IoActor.MaxAttemptsNumber) None else Option(backoff.backoffMillis milliseconds)

  /**
    * Json batch call back for a batched request
    */
  lazy val callback: JsonBatchCallback[U] = new JsonBatchCallback[U]() {
    def onSuccess(response: U, httpHeaders: HttpHeaders) = onSuccessCallback(response, httpHeaders)
    def onFailure(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders) = onFailureCallback(googleJsonError, httpHeaders)
  }

  /**
    * Increment backoff time and attempt count
    */
  lazy val next: GcsBatchCommandContext[T, U] = {
    this.copy(backoff = backoff.next, currentAttempt = currentAttempt + 1, promise = Promise[BatchResponse])
  }

  /**
    * Only increment backoff. To be used for failure thas should be retried infinitely
    */
  lazy val nextTransient: GcsBatchCommandContext[T, U] = {
    this.copy(backoff = backoff.next, promise = Promise[BatchResponse])
  }

  /**
    * Queue the request for batching
    */
  def queue(batchRequest: BatchRequest) = request.operation.queue(batchRequest, callback)

  /**
    * On success callback. Transform the request response to a stream-ready response that can complete the promise
    */
  private def onSuccessCallback(response: U, httpHeaders: HttpHeaders) = handleSuccessOrNextRequest(request.onSuccess(response, httpHeaders))

  private def handleSuccessOrNextRequest(successResult: Either[T, GcsBatchIoCommand[T, U]]) = {
    val promiseResponse: BatchResponse = successResult match {
      // Left means the command is complete, so just create the corresponding IoSuccess with the value
      case Left(responseValue) => Left(success(responseValue))
      // Right means there is a subsequent request to be executed, clone this context with the new request and a new promise
      case Right(nextCommand) => Right(this.copy(request = nextCommand, promise = Promise[BatchResponse]))
    }

    promise.trySuccess(promiseResponse)
    ()
  }

  /**
    * On failure callback. Fail the promise with a StorageException
    */
  private def onFailureCallback(googleJsonError: GoogleJsonError, httpHeaders: HttpHeaders) = {
    if (isProjectNotProvidedError(googleJsonError)) {
      // Returning an Either.Right here means that the operation is not complete and that we need to do another request
      handleSuccessOrNextRequest(Right(request.withUserProject))
    } else {
      (request.onFailure(googleJsonError, httpHeaders), request) match {
        case (Some(successValue), _) => handleSuccessOrNextRequest(successValue)
        case (None, singleFile: SingleFileIoCommand[_]) if googleJsonError.getCode == 404 =>
          // Make the message clearer if it's a file not found error
          googleJsonError.setMessage(s"Object ${singleFile.file.pathAsString} does not exist")
          promise.tryFailure(new FileNotFoundException(singleFile.file.pathAsString))
        case (None, _) => promise.tryFailure(new StorageException(googleJsonError))
      }
      ()
    }
  }
}
