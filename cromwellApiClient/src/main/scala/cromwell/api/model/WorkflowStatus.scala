package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

// ********* !!!!!!!!!! ********
//
// WARNING! This is a Cromwell API class. If you aren't changing the API client, you probably
// want to look elsewhere (maybe in cromwell.core?)
//
// ********* !!!!!!!!!! ********

/**
  * ADT tree to describe Cromwell workflow statuses, both terminal and non-terminal
  */
sealed trait WorkflowStatus

sealed trait TerminalStatus extends WorkflowStatus
case object Aborted extends TerminalStatus
case object Failed extends TerminalStatus
case object Succeeded extends TerminalStatus

sealed trait NonTerminalStatus extends WorkflowStatus
case object Submitted extends NonTerminalStatus
case object Running extends NonTerminalStatus
case object Aborting extends NonTerminalStatus

object WorkflowStatus {
  def apply(status: String): WorkflowStatus = {
    status match {
      case "Submitted" => Submitted
      case "Running" => Running
      case "Aborting" => Aborting
      case "Aborted" => Aborted
      case "Failed" => Failed
      case "Succeeded" => Succeeded
      case bad => throw new IllegalArgumentException(s"No such status: $bad")
    }
  }

  def apply(workflowStatus: CromwellStatus): WorkflowStatus = apply(workflowStatus.status)
}

object WorkflowStatusJsonFormatter extends DefaultJsonProtocol {
  implicit object WorkflowStatusJsonFormat extends RootJsonFormat[WorkflowStatus] {
    def write(status: WorkflowStatus) = new JsString(status.toString)
    def read(value: JsValue) = value match {
      case JsString(string) => WorkflowStatus(string)
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into a WorkflowStatus")
    }
  }
}
