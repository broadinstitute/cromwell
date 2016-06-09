package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, Starting, Running, Failed, Preempted, Done, Aborted = Value

  implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
    def isTerminal: Boolean = {
      Set(Failed, Done, Aborted, Preempted) contains status
    }
  }

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
  }
}
