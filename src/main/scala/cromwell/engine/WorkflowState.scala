package cromwell.engine

sealed trait WorkflowState {
  def isTerminal: Boolean
}

object WorkflowState {
  private lazy val WorkflowState = Seq(WorkflowSubmitted, WorkflowRunning, WorkflowFailed, WorkflowSucceeded, WorkflowAborting, WorkflowAborted)

  def fromString(str: String): WorkflowState = WorkflowState.find(_.toString == str).getOrElse(
    throw new NoSuchElementException(s"No such WorkflowState: $str"))
}

case object WorkflowSubmitted extends WorkflowState {
  override val toString: String = "Submitted"
  override val isTerminal = false
}

case object WorkflowRunning extends WorkflowState {
  override val toString: String = "Running"
  override val isTerminal = false
}

case object WorkflowAborting extends WorkflowState {
  override val toString: String = "Aborting"
  override val isTerminal = false
}

case object WorkflowFailed extends WorkflowState {
  override val toString: String = "Failed"
  override val isTerminal = true
}

case object WorkflowSucceeded extends WorkflowState {
  override val toString: String = "Succeeded"
  override val isTerminal = true
}

case object WorkflowAborted extends WorkflowState {
  override val toString: String = "Aborted"
  override val isTerminal = true
}