package cromwell.backend.google.batch.monitoring

import cats.data.NonEmptyList
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.instrumentation.CromwellInstrumentation._

object BatchInstrumentation {
  private val BatchKey = NonEmptyList.of("batch")
  private val BatchPollKey = BatchKey.concatNel("poll")
  private val BatchRunKey = BatchKey.concatNel("run")
  private val BatchAbortKey = BatchKey.concatNel("abort")

  private val BatchPollFailedKey = BatchPollKey.concatNel(FailureKey)
  private val BatchRunFailedKey = BatchRunKey.concatNel(FailureKey)
  private val BatchAbortFailedKey = BatchAbortKey.concatNel(FailureKey)
  implicit class StatsDPathGoogleEnhanced(val statsDPath: InstrumentationPath) extends AnyVal {
    def withGoogleThrowable(failure: Throwable) =
      statsDPath.withThrowable(failure, statusCodeExtractor)
  }

  private def statusCodeExtractor(ex: Throwable): Option[Int] = ex match {
    case batchEx: BatchApiException =>
      batchEx.cause match {
        case GoogleBatchException(apiException) => Option(apiException.getStatusCode.getCode.getHttpStatusCode)
        case _ => None
      }
    case _ => None
  }
}

trait BatchInstrumentation extends CromwellInstrumentation {
  import BatchInstrumentation._

  def pollSuccess(): Unit = increment(BatchPollKey.concatNel(SuccessKey), BackendPrefix)
  def runSuccess(): Unit = increment(BatchRunKey.concatNel(SuccessKey), BackendPrefix)
  def abortSuccess(): Unit = increment(BatchAbortKey.concatNel(SuccessKey), BackendPrefix)

  def failedQuery(failedQuery: BatchApiRequestFailed): Unit = failedQuery.query match {
    case _: BatchStatusPollRequest =>
      increment(BatchPollFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: BatchRunCreationRequest =>
      increment(BatchRunFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: BatchAbortRequest =>
      increment(BatchAbortFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
  }

  def retriedQuery(failedQuery: BatchApiRequestFailed): Unit = failedQuery.query match {
    case _: BatchStatusPollRequest =>
      increment(BatchPollFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: BatchRunCreationRequest =>
      increment(BatchRunFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: BatchAbortRequest =>
      increment(BatchAbortFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
  }

  def updateQueueSize(size: Int) = sendGauge(BatchKey.concatNel("queue_size"), size.toLong, BackendPrefix)

}
