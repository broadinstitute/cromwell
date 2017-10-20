package cromwell.api.model

import java.util.UUID

import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

// ********* !!!!!!!!!! ********
//
// WARNING! This is the Cromwell API version of WorkflowId. If you aren't changing the API client, you probably
// want cromwell.core.WorkflowId instead!
//
// ********* !!!!!!!!!! ********

final case class WorkflowId(id: UUID) extends AnyVal {
  override def toString = id.toString
  def shortString = id.toString.split("-")(0)
}

object WorkflowId {
  def fromString(id: String): WorkflowId = new WorkflowId(UUID.fromString(id))
  def randomId() = WorkflowId(UUID.randomUUID())
}

object WorkflowIdJsonFormatter extends DefaultJsonProtocol {
  implicit object WorkflowIdJsonFormat extends RootJsonFormat[WorkflowId] {
    def write(id: WorkflowId) = JsString(id.id.toString)
    def read(value: JsValue) = value match {
      case JsString(s) => WorkflowId.fromString(s)
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into a ShardIndex")
    }
  }
}

