package cromwell.backend.google.pipelines.common

import cats.data.NonEmptyList
import cromwell.backend.google.pipelines.common.PapiInstrumentation._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.filesystems.gcs.GoogleUtil
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.instrumentation.CromwellInstrumentation._

object PapiInstrumentation {
  private val PapiKey = NonEmptyList.of("papi")
  private val PapiPollKey = PapiKey.concatNel("poll")
  private val PapiRunKey = PapiKey.concatNel("run")
  private val PapiAbortKey = PapiKey.concatNel("abort")

  private val PapiPollFailedKey = PapiPollKey.concatNel(FailureKey)
  private val PapiRunFailedKey = PapiRunKey.concatNel(FailureKey)
  private val PapiAbortFailedKey = PapiAbortKey.concatNel(FailureKey)
  private val PapiPollRetriedKey = PapiPollKey.concatNel(RetryKey)
  private val PapiRunRetriedKey = PapiRunKey.concatNel(RetryKey)
  private val PapiAbortRetriedKey = PapiAbortKey.concatNel(RetryKey)

  implicit class StatsDPathGoogleEnhanced(val statsDPath: InstrumentationPath) extends AnyVal {
    def withGoogleThrowable(failure: Throwable) = {
      statsDPath.withThrowable(failure, GoogleUtil.extractStatusCode)
    }
  }
}

trait PapiInstrumentation extends CromwellInstrumentation {
  def pollSuccess() = increment(PapiPollKey.concatNel(SuccessKey), BackendPrefix)
  def runSuccess() = increment(PapiRunKey.concatNel(SuccessKey), BackendPrefix)
  def abortSuccess() = increment(PapiAbortKey.concatNel(SuccessKey), BackendPrefix)

  def failedQuery(failedQuery: PAPIApiRequestFailed) = failedQuery.query match {
    case _: PAPIStatusPollRequest => increment(PapiPollFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: PAPIRunCreationRequest => increment(PapiRunFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: PAPIAbortRequest => increment(PapiAbortFailedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
  }

  def retriedQuery(failedQuery: PAPIApiRequestFailed) = failedQuery.query match {
    case _: PAPIStatusPollRequest => increment(PapiPollRetriedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: PAPIRunCreationRequest => increment(PapiRunRetriedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
    case _: PAPIAbortRequest => increment(PapiAbortRetriedKey.withGoogleThrowable(failedQuery.cause.cause), BackendPrefix)
  }

  def updateQueueSize(size: Int) = sendGauge(PapiKey.concatNel("queue_size"), size.toLong, BackendPrefix)
}
