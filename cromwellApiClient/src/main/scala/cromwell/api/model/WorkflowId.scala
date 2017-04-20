package cromwell.api.model

import java.util.UUID

// ********* !!!!!!!!!! ********
//
// WARNING! This is the Cromwell API version of WorkflowId. If you aren't changing the API client, you probably
// want cromwell.core.WorkflowId instead!
//
// ********* !!!!!!!!!! ********

case class WorkflowId(id: UUID) {
  override def toString = id.toString
  def shortString = id.toString.split("-")(0)
}

object WorkflowId {
  def fromString(id: String): WorkflowId = new WorkflowId(UUID.fromString(id))
  def randomId() = WorkflowId(UUID.randomUUID())
}
