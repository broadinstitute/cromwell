package cromwell.api.model

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

object CromwellBackendsJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellBackendsFormat: RootJsonFormat[CromwellBackends] = jsonFormat2(CromwellBackends)
}

final case class CromwellBackends(defaultBackend: String, supportedBackends: List[String])
