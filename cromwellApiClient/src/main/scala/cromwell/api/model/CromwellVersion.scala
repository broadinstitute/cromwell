package cromwell.api.model

import spray.json.DefaultJsonProtocol

object CromwellVersionJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellVersionFormat = jsonFormat1(CromwellVersion)
}

case class CromwellVersion(cromwell: String)
