package cromwell.core

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, QueuedInCromwell, Starting, Running, Aborting, Failed, RetryableFailure, Done, Bypassed, Aborted, Unstartable = Value
  val TerminalStatuses = Set(Failed, Done, Aborted, Bypassed, Unstartable)
  val TerminalOrRetryableStatuses = TerminalStatuses + RetryableFailure
  val NonTerminalStatuses = values.diff(TerminalOrRetryableStatuses)
  val ActiveStatuses = Set(QueuedInCromwell, Starting, Running, Aborting)

  implicit val ExecutionStatusOrdering = Ordering.by { status: ExecutionStatus =>
    status match {
      case NotStarted => 0
      case QueuedInCromwell => 1
      case Starting => 2
      case Running => 3
      case Aborting => 4
      case Unstartable => 5
      case Aborted => 6
      case Bypassed => 7
      case RetryableFailure => 8
      case Failed => 9
      case Done => 10
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
