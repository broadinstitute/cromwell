package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, WaitingForQueueSpace, QueuedInCromwell, Starting, Running, Aborting, Failed, RetryableFailure, Done, Bypassed, Aborted, Unstartable = Value
  val TerminalStatuses = Set(Failed, Done, Aborted, Bypassed, Unstartable)
  val TerminalOrRetryableStatuses = TerminalStatuses + RetryableFailure
  val NonTerminalStatuses = values.diff(TerminalOrRetryableStatuses)
  val ActiveStatuses = Set(WaitingForQueueSpace, QueuedInCromwell, Starting, Running, Aborting)

  implicit val ExecutionStatusOrdering = Ordering.by { status: ExecutionStatus =>
    status match {
      case NotStarted => 0
      case WaitingForQueueSpace => 1
      case QueuedInCromwell => 2
      case Starting => 3
      case Running => 4
      case Aborting => 5
      case Unstartable => 6
      case Aborted => 7
      case Bypassed => 8
      case RetryableFailure => 9
      case Failed => 10
      case Done => 11
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
