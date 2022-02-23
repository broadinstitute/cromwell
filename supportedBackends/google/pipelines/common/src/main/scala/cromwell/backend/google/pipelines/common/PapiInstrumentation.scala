package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.PapiInstrumentation._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.filesystems.gcs.GoogleUtil
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.instrumentation.CromwellInstrumentation._

object PapiInstrumentation {
  private val PapiKey = InstrumentationPath.withParts("papi")
  private val PapiPollKey = PapiKey.withParts("poll")
  private val PapiRunKey = PapiKey.withParts("run")
  private val PapiAbortKey = PapiKey.withParts("abort")

  private val PapiPollFailedKey = PapiPollKey.withParts(FailureKey)
  private val PapiRunFailedKey = PapiRunKey.withParts(FailureKey)
  private val PapiAbortFailedKey = PapiAbortKey.withParts(FailureKey)
  private val PapiPollRetriedKey = PapiPollKey.withParts(RetryKey)
  private val PapiRunRetriedKey = PapiRunKey.withParts(RetryKey)
  private val PapiAbortRetriedKey = PapiAbortKey.withParts(RetryKey)

  def pathFromFailedQuery(path: InstrumentationPath, failedQuery: PAPIApiRequestFailed): InstrumentationPath =
    path.withThrowable(failedQuery.cause.cause, GoogleUtil.extractStatusCode)
}

trait PapiInstrumentation extends CromwellInstrumentation {
  def pollSuccess() = increment(PapiPollKey.withParts(SuccessKey), BackendPrefix)
  def runSuccess() = increment(PapiRunKey.withParts(SuccessKey), BackendPrefix)
  def abortSuccess() = increment(PapiAbortKey.withParts(SuccessKey), BackendPrefix)

  def failedQuery(failedQuery: PAPIApiRequestFailed) = failedQuery.query match {
    case _: PAPIStatusPollRequest => increment(pathFromFailedQuery(PapiPollFailedKey, failedQuery), BackendPrefix)
    case _: PAPIRunCreationRequest => increment(pathFromFailedQuery(PapiRunFailedKey, failedQuery), BackendPrefix)
    case _: PAPIAbortRequest => increment(pathFromFailedQuery(PapiAbortFailedKey, failedQuery), BackendPrefix)
  }

  def retriedQuery(failedQuery: PAPIApiRequestFailed) = failedQuery.query match {
    case _: PAPIStatusPollRequest => increment(pathFromFailedQuery(PapiPollRetriedKey, failedQuery), BackendPrefix)
    case _: PAPIRunCreationRequest => increment(pathFromFailedQuery(PapiRunRetriedKey, failedQuery), BackendPrefix)
    case _: PAPIAbortRequest => increment(pathFromFailedQuery(PapiAbortRetriedKey, failedQuery), BackendPrefix)
  }

  def updateQueueSize(size: Int) = sendGauge(PapiKey.withParts("queue_size"), size.toLong, BackendPrefix)
}
