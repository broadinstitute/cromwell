package cromwell.backend.google.pipelines.batch.monitoring

import cats.data.NonEmptyList
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
}

// TODO: We could propagate details about the failed requests like lifesciences backend does (e.g. response code)
trait BatchInstrumentation extends CromwellInstrumentation {
  import BatchInstrumentation._

  def pollSuccess() = increment(BatchPollKey.concatNel(SuccessKey), BackendPrefix)
  def pollFailed() = increment(BatchPollFailedKey.concatNel(FailureKey), BackendPrefix)

  def runSuccess() = increment(BatchRunKey.concatNel(SuccessKey), BackendPrefix)
  def runFailed() = increment(BatchRunFailedKey.concatNel(FailureKey), BackendPrefix)

  def abortSuccess() = increment(BatchAbortKey.concatNel(SuccessKey), BackendPrefix)
  def abortFailed() = increment(BatchAbortFailedKey.concatNel(FailureKey), BackendPrefix)

}
