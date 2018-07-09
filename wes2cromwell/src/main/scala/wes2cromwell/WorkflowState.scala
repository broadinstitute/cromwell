package wes2cromwell

import spray.json.{ JsString, JsValue, RootJsonFormat }

sealed trait WorkflowState

object WorkflowState {
  lazy val WorkflowStateValues = Seq(UNKNOWN, QUEUED, INITIALIZING, RUNNING, PAUSED, COMPLETE, EXECUTOR_ERROR, SYSTEM_ERROR, CANCELED)

  implicit object WorkflowStateFormat extends RootJsonFormat[WorkflowState] {
    def write(obj: WorkflowState): JsValue = JsString(obj.toString)

    // TODO: not sure this works and nothing in the code tests it right now
    def read(json: JsValue): WorkflowState = {
      json match {
        case x: JsString => WorkflowState.withName(x.value)
        case _ => throw new IllegalArgumentException("Expected a string value for state")
      }
    }
  }

  def withName(str: String): WorkflowState = WorkflowStateValues.find(_.toString == str).getOrElse(
    throw new NoSuchElementException(s"No such WorkflowState: $str")
  )

  case object UNKNOWN extends WorkflowState {
    override val toString: String = "UNKNOWN"
  }
  case object QUEUED extends WorkflowState {
    override val toString: String = "QUEUED"
  }
  case object INITIALIZING extends WorkflowState {
    override val toString: String = "INITIALIZING"
  }
  case object RUNNING extends WorkflowState {
    override val toString = "RUNNING"
  }
  case object PAUSED extends WorkflowState {
    override val toString = "PAUSED"
  }
  case object COMPLETE extends WorkflowState {
    override val toString = "COMPLETE"
  }
  case object EXECUTOR_ERROR extends WorkflowState {
    override val toString = "EXECUTOR_ERROR"
  }
  case object SYSTEM_ERROR extends WorkflowState {
    override val toString = "SYSTEM_ERROR"
  }
  case object CANCELED extends WorkflowState {
    override val toString = "CANCELED"
  }
}

