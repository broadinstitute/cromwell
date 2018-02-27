package cromwell.backend.impl.jes.statuspolling

import akka.actor.Actor
import cats.data.NonEmptyList
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiQueryFailed, JesRunCreationQuery, JesStatusPollQuery}
import cromwell.backend.impl.jes.statuspolling.PapiInstrumentation._
import cromwell.backend.instrumentation.BackendInstrumentation._
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.filesystems.gcs.GoogleUtil
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentationActor

object PapiInstrumentation {
  private val PapiKey = NonEmptyList.of("papi")
  private val PapiPollKey = PapiKey.concat("poll")
  private val PapiRunKey = PapiKey.concat("run")
  
  private val PapiPollFailedKey = PapiPollKey.concat(FailureKey)
  private val PapiRunFailedKey = PapiRunKey.concat(FailureKey)
  private val PapiPollRetriedKey = PapiPollKey.concat(RetryKey)
  private val PapiRunRetriedKey = PapiRunKey.concat(RetryKey)
  
  implicit class StatsDPathGoogleEnhanced(val statsDPath: InstrumentationPath) extends AnyVal {
    def withGoogleThrowable(failure: Throwable) = {
      statsDPath.withThrowable(failure, GoogleUtil.extractStatusCode)
    }
  }
}

trait PapiInstrumentation extends CromwellInstrumentationActor { this: Actor =>
  def pollSuccess() = increment(PapiPollKey.concat(SuccessKey))
  def runSuccess() = increment(PapiRunKey.concat(SuccessKey))

  def failedQuery(failedQuery: JesApiQueryFailed) = failedQuery.query match {
    case _: JesStatusPollQuery => increment(PapiPollFailedKey.withGoogleThrowable(failedQuery.cause.e), BackendPrefix)
    case _: JesRunCreationQuery => increment(PapiRunFailedKey.withGoogleThrowable(failedQuery.cause.e), BackendPrefix)
  }

  def retriedQuery(failedQuery: JesApiQueryFailed) = failedQuery.query match {
    case _: JesStatusPollQuery => increment(PapiPollRetriedKey.withGoogleThrowable(failedQuery.cause.e), BackendPrefix)
    case _: JesRunCreationQuery => increment(PapiRunRetriedKey.withGoogleThrowable(failedQuery.cause.e), BackendPrefix)
  }
  
  def updateQueueSize(size: Int) = sendGauge(PapiKey.concat("queue_size"), size.toLong, BackendPrefix)
}
