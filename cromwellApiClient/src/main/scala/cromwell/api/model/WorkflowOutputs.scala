package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsValue}
import spray.json.RootJsonFormat

object WorkflowOutputsJsonSupport extends DefaultJsonProtocol {
  implicit val OutputResponseFormat: RootJsonFormat[WorkflowOutputs] = jsonFormat2(WorkflowOutputs)
}

case class WorkflowOutputs(id: String, outputs: JsValue)
