package centaur

import spray.json.DefaultJsonProtocol

object CromwellStatusJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellStatusFormat = jsonFormat2(CromwellStatus)
}

case class CromwellStatus(id: String, status: String)