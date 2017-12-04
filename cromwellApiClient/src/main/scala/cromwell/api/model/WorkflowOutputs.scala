package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsValue}

object WorkflowOutputsJsonSupport extends DefaultJsonProtocol {
  implicit val OutputResponseFormat = jsonFormat2(WorkflowOutputs)
}

case class WorkflowOutputs(id: String, outputs: JsValue)
