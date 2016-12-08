package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, QueuedInCromwell, Starting, Running, Failed, Preempted, Done, Bypassed, Aborted = Value
  val TerminalStatuses = Set(Failed, Done, Aborted, Preempted, Bypassed)

  implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
    def isTerminal: Boolean = {
      TerminalStatuses contains status
    }

    def isDoneOrBypassed: Boolean = status == Done || status == Bypassed
  }

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
  }
}
