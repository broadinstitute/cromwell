package cromwell.engine.instrumentation

import akka.actor.Actor
import cromwell.backend.BackendJobExecutionActor._
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.engine.instrumentation.JobInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.duration.FiniteDuration

object JobInstrumentation {
  private def backendJobExecutionResponsePaths(response: BackendJobExecutionResponse) = response match {
    case _: JobSucceededResponse => SuccessKey
    case _: JobAbortedResponse => AbortedKey
    case _: JobFailedNonRetryableResponse => FailureKey
    case _: JobFailedRetryableResponse => RetryKey
  }
}

/**
  * Provides helper methods for Job instrumentation
  */
trait JobInstrumentation extends CromwellInstrumentationActor { this: Actor =>

  /**
    * Generic method to increment a workflow related counter metric value
    */
  def incrementJob(statsDPath: InstrumentationPath) = increment(statsDPath, JobPrefix)

  /**
    * Generic method to add a workflow related timing metric value
    */
  def setTimingJob(statsDPath: InstrumentationPath, duration: FiniteDuration): Unit = {
    sendTiming(statsDPath, duration, JobPrefix)
  }

  /**
    * Generic method to update a job related gauge metric value
    */
  def sendGaugeJob(statsDPath: InstrumentationPath, value: Long): Unit = {
    sendGauge(statsDPath, value, JobPrefix)
  }

  /**
    * Add a timing value for the run time of a job in a given state
    */
  def setJobTimePerState(response: BackendJobExecutionResponse, duration: FiniteDuration): Unit = {
    setTimingJob(backendJobExecutionResponsePaths(response), duration)
  }
}
