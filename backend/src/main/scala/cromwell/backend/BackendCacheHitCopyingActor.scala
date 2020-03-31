package cromwell.backend

import cromwell.core.JobKey
import cromwell.core.simpleton.WomValueSimpleton

object BackendCacheHitCopyingActor {
  final case class CopyOutputsCommand(womValueSimpletons: Seq[WomValueSimpleton], jobDetritusFiles: Map[String, String], returnCode: Option[Int])

  final case class CopyingOutputsFailedResponse(jobKey: JobKey, cacheCopyAttempt: Int, failure: LoggableOrMetricableCacheCopyError)

  sealed trait LoggableOrMetricableCacheCopyError
  final case class LoggableCacheCopyError(failure: Throwable) extends LoggableOrMetricableCacheCopyError
  final case class MetricableCacheCopyError(failureCategory: String) extends LoggableOrMetricableCacheCopyError
}
