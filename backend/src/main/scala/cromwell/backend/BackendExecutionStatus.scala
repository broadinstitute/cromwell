package cromwell.backend

object BackendExecutionStatus extends Enumeration {
  type BackendExecutionStatus = Value
  val NotStarted, Starting, Running, Failed, Preempted, Done, Aborted = Value
}
