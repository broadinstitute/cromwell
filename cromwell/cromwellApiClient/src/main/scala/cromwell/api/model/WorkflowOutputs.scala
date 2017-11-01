package cromwell.api.model

import spray.json.DefaultJsonProtocol

object WorkflowOutputsJsonSupport extends DefaultJsonProtocol {
  implicit val OutputResponseFormat = jsonFormat2(WorkflowOutputs)
}

case class WorkflowOutputs(id: String, outputs: Map[String, String])
