package cromwell.backend.google.pipelines.batch.runnable

sealed trait RunnableFlag

object RunnableFlag {
//  case object Unspecified extends RunnableFlag
  case object IgnoreExitStatus extends RunnableFlag
  case object RunInBackground extends RunnableFlag
  case object AlwaysRun extends RunnableFlag
//  case object EnableFuse extends RunnableFlag
//  case object PublishExposedPorts extends RunnableFlag
//  case object DisableImagePrefetch extends RunnableFlag
}
