package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, QueuedInCromwell, Starting, Running, Failed, Preempted, Done, Aborted = Value
  val TerminalStatuses = Set(Failed, Done, Aborted, Preempted)

  implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
    def isTerminal: Boolean = {
      TerminalStatuses contains status
    }
  }

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
  }
}
