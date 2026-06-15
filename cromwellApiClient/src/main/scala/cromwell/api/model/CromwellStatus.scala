package cromwell.api.model

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

object CromwellStatusJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellStatusFormat: RootJsonFormat[CromwellStatus] = jsonFormat2(CromwellStatus)
}

case class CromwellStatus(id: String, status: String)
