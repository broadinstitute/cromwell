package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, QueuedInCromwell, Starting, Running, Failed, RetryableFailure, Done, Bypassed, Aborted = Value
  val TerminalStatuses = Set(Failed, Done, Aborted, Bypassed)
  val TerminalOrRetryableStatuses = TerminalStatuses + RetryableFailure

  implicit val ExecutionStatusOrdering = Ordering.by { status: ExecutionStatus =>
    status match {
      case NotStarted => 0
      case QueuedInCromwell => 1
      case Starting => 2
      case Running => 3
      case Aborted => 4
      case Bypassed => 5
      case RetryableFailure => 6
      case Failed => 7
      case Done => 8
    }
  }
  
  implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
    def isTerminal: Boolean = TerminalStatuses contains status

    def isTerminalOrRetryable: Boolean = TerminalOrRetryableStatuses contains status

    def isDoneOrBypassed: Boolean = status == Done || status == Bypassed
  }

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
  }
}
