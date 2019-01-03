package cromwell.core

import cats.Semigroup


sealed trait WorkflowState {
  def isTerminal: Boolean
  protected def ordinal: Int
  def combine(that: WorkflowState): WorkflowState = if (this.ordinal > that.ordinal) this else that
}

object WorkflowState {
  lazy val WorkflowStateValues = Seq(WorkflowOnHold, WorkflowSubmitted, WorkflowRunning, WorkflowFailed, WorkflowSucceeded, WorkflowAborting, WorkflowAborted)

  def withName(str: String): WorkflowState = WorkflowStateValues.find(_.toString.equalsIgnoreCase(str)).getOrElse(
    throw new NoSuchElementException(s"No such WorkflowState: $str"))

  implicit val WorkflowStateSemigroup = new Semigroup[WorkflowState] {
    override def combine(f1: WorkflowState, f2: WorkflowState): WorkflowState = f1.combine(f2)
  }

  implicit val WorkflowStateOrdering = Ordering.by { self: WorkflowState => self.ordinal }
}

case object WorkflowOnHold extends WorkflowState{
  override val toString: String = "On Hold"
  override val isTerminal = false
  override val ordinal = 0
}

case object WorkflowSubmitted extends WorkflowState {
  override val toString: String = "Submitted"
  override val isTerminal = false
  override val ordinal = 1
}

case object WorkflowRunning extends WorkflowState {
  override val toString: String = "Running"
  override val isTerminal = false
  override val ordinal = 2
}

case object WorkflowAborting extends WorkflowState {
  override val toString: String = "Aborting"
  override val isTerminal = false
  override val ordinal = 3
}

case object WorkflowAborted extends WorkflowState {
  override val toString: String = "Aborted"
  override val isTerminal = true
  override val ordinal = 4
}

case object WorkflowSucceeded extends WorkflowState {
  override val toString: String = "Succeeded"
  override val isTerminal = true
  override val ordinal = 5
}

case object WorkflowFailed extends WorkflowState {
  override val toString: String = "Failed"
  override val isTerminal = true
  override val ordinal = 6
}
