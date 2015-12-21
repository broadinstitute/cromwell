package cromwell.engine

object ExecutionStatus extends Enumeration {
  type ExecutionStatus = Value
  val NotStarted, Starting, Running, Failed, Done, Aborted = Value

  implicit class EnhancedExecutionStatus(val status: ExecutionStatus) extends AnyVal {
    def isTerminal: Boolean = {
      Seq(Failed, Done, Aborted) contains status
    }
  }

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toExecutionStatus: ExecutionStatus = ExecutionStatus.withName(string)
  }
}