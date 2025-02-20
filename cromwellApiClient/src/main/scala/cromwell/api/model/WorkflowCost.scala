package cromwell.api.model

import spray.json.DefaultJsonProtocol

object WorkflowCostJsonSupport extends DefaultJsonProtocol {
  implicit val WorkflowCostJsonFormat = jsonFormat5(WorkflowCost)
}

final case class WorkflowCost(errors: List[String], id: String, cost: BigDecimal, status: String, currency: String)
