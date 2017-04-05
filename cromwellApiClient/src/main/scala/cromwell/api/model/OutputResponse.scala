package cromwell.api.model

import spray.json.DefaultJsonProtocol

object OutputResponseJsonSupport extends DefaultJsonProtocol {
  implicit val OutputResponseFormat = jsonFormat2(OutputResponse)
}

case class OutputResponse(id: String, outputs: Map[String, String])
