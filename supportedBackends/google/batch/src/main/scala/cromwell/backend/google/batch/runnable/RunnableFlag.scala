package cromwell.backend.google.batch.runnable

sealed trait RunnableFlag

object RunnableFlag {
  case object IgnoreExitStatus extends RunnableFlag
  case object RunInBackground extends RunnableFlag
  case object AlwaysRun extends RunnableFlag
}
