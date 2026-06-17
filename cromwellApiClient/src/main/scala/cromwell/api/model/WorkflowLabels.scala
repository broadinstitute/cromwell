package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsObject}
import spray.json.RootJsonFormat

object WorkflowLabelsJsonSupport extends DefaultJsonProtocol {
  implicit val LabelsResponseFormat: RootJsonFormat[WorkflowLabels] = jsonFormat2(WorkflowLabels)
}

final case class WorkflowLabels(id: String, labels: JsObject)
