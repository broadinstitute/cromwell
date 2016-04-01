package cromwell.core

import java.util.UUID

case class WorkflowId(id: UUID) {
  override def toString = id.toString
  def shortString = id.toString.split("-")(0)
}

object WorkflowId {
  def fromString(id: String): WorkflowId = new WorkflowId(UUID.fromString(id))
  def randomId() = WorkflowId(UUID.randomUUID())
}