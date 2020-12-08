package cromwell.engine.io

import java.io.IOException
import java.net.{SocketException, SocketTimeoutException}

import com.google.cloud.storage.StorageException
import cromwell.engine.io.gcs.GcsBatchFlow.BatchFailedException
import javax.net.ssl.SSLException

object RetryableRequestSupport {

  def isFatal(failure: Throwable): Boolean = !isRetryable(failure)

  /**
    * Failures that are considered retryable.
    * Retrying them should increase the "retry counter"
    */
  def isRetryable(failure: Throwable): Boolean = failure match {
    case gcs: StorageException => gcs.isRetryable ||
      isRetryable(gcs.getCause) ||
      AdditionalRetryableHttpCodes.contains(gcs.getCode) ||
      Option(gcs.getMessage).exists(msg =>
        AdditionalRetryableErrorMessages.contains(msg.toLowerCase))
    case _: SSLException => true
    case _: BatchFailedException => true
    case _: SocketException => true
    case _: SocketTimeoutException => true
    case ioE: IOException if Option(ioE.getMessage).exists(_.contains("Error getting access token for service account")) => true
    case ioE: IOException => isGcs500(ioE) || isGcs503(ioE) || isGcs504(ioE)
    case other => isTransient(other)
  }

  private val AdditionalRetryableHttpCodes = List(
    // HTTP 410: Gone
    // From Google doc (https://cloud.google.com/storage/docs/json_api/v1/status-codes):
    // "You have attempted to use a resumable upload session that is no longer available.
    // If the reported status code was not successful and you still wish to upload the file, you must start a new session."
    410,
    // Some 503 errors seem to yield "false" on the "isRetryable" method because they are not retried.
    // The CloudStorage exception mechanism is not flawless yet (https://github.com/GoogleCloudPlatform/google-cloud-java/issues/1545)
    // so that could be the cause.
    // For now explicitly lists 503 as a retryable code here to work around that.
    503
  )

  // Error messages not included in the list of built-in GCS retryable errors (com.google.cloud.storage.StorageException) but that we still want to retry
  private val AdditionalRetryableErrorMessages = List(
    "Connection closed prematurely"
  ).map(_.toLowerCase)

  /**
    * ATTENTION: Transient failures are retried *forever*
    * Be careful when adding error codes to this method.
    * Currently only 429 (= quota exceeded are considered truly transient)
    */
  def isTransient(failure: Throwable): Boolean = failure match {
    case gcs: StorageException => gcs.getCode == 429
    case _ => false
  }

  def isGcs500(failure: Throwable): Boolean = {
    val serverErrorPattern = ".*Could not read from gs.+500 Internal Server Error.*"
    Option(failure.getMessage).exists(_.matches(serverErrorPattern))
  }

  def isGcs503(failure: Throwable): Boolean = {
    val serverErrorPattern = ".*Could not read from gs.+503 Service Unavailable.*"
    Option(failure.getMessage).exists(_.matches(serverErrorPattern))
  }

  def isGcs504(failure: Throwable): Boolean = {
    val serverErrorPattern = ".*Could not read from gs.+504 Gateway Timeout.*"
    Option(failure.getMessage).exists(_.matches(serverErrorPattern))
  }
}
