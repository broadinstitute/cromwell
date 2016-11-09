package centaur.api
import spray.json.DefaultJsonProtocol

import scala.util.{Failure, Success, Try}

object CromwellBackendsJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellBackendsFormat = jsonFormat2(CromwellBackends)
}

final case class CromwellBackends(defaultBackend: String, supportedBackends: List[String])

// Stupid Spray
object CromwellBackendsCompanion {
  def supportedBackends: List[String] = {
    CromwellClient.backends match {
      case Success(b) => b.supportedBackends map { _.toLowerCase }
      case Failure(e) => throw e
    }
  }
}
