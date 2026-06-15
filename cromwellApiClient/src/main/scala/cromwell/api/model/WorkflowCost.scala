package cromwell.api.model

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

object WorkflowCostJsonSupport extends DefaultJsonProtocol {
  implicit val WorkflowCostJsonFormat: RootJsonFormat[WorkflowCost] = jsonFormat5(WorkflowCost)
}

final case class WorkflowCost(errors: List[String], id: String, cost: BigDecimal, status: String, currency: String)
