package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsObject}

object WorkflowLabelsJsonSupport extends DefaultJsonProtocol {
  implicit val LabelsResponseFormat = jsonFormat2(WorkflowLabels)
}

final case class WorkflowLabels(id: String, labels: JsObject)
