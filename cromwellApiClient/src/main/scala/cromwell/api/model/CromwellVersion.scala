package cromwell.api.model

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

object CromwellVersionJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellVersionFormat: RootJsonFormat[CromwellVersion] = jsonFormat1(CromwellVersion)
}

case class CromwellVersion(cromwell: String)
