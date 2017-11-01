package cromwell.api.model

import spray.json.DefaultJsonProtocol

object CromwellBackendsJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellBackendsFormat = jsonFormat2(CromwellBackends)
}

final case class CromwellBackends(defaultBackend: String, supportedBackends: List[String])
