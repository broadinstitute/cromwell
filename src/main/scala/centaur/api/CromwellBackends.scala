package centaur.api
import centaur.CentaurConfig
import CromwellClient._
import centaur.api.CromwellBackendsJsonSupport.CromwellBackendsFormat
import spray.client.pipelining._
import spray.json.DefaultJsonProtocol
import spray.httpx.SprayJsonSupport._

import scala.util.{Failure, Success, Try}


object CromwellBackendsJsonSupport extends DefaultJsonProtocol {
  implicit val CromwellBackendsFormat = jsonFormat2(CromwellBackends)
}

final case class CromwellBackends(defaultBackend: String, supportedBackends: List[String])

// Stupid Spray
object CromwellBackendsCompanion {
  def cromwellBackends: Try[CromwellBackends] = {
    val response = Pipeline[CromwellBackends].apply(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/backends"))
    sendReceiveFutureCompletion(response)
  }

  def supportedBackends: List[String] = {
    cromwellBackends match {
      case Success(b) => b.supportedBackends map { _.toLowerCase }
      case Failure(e) => throw e
    }
  }
}
