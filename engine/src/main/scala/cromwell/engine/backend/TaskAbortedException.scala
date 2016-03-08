package cromwell.engine.backend

class TaskAbortedException extends Exception {
  override def getMessage: String = "The task was aborted."
}
