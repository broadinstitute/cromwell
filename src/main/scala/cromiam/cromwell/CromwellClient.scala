package cromiam.cromwell

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import com.softwaremill.sttp._
import cromiam.cromwell.CromwellClient._
import cromiam.server.status.StatusCheckedSubsystem
import cromwell.api.model.CromwellStatus
import cromwell.api.model.CromwellStatusJsonSupport._
import spray.json._

import scala.concurrent.{ExecutionContextExecutor, Future}

/**
  * Provides a CromIAM specific handle for Cromwell communication
  */
class CromwellClient(scheme: String, interface: String, port: Int)(implicit system: ActorSystem,
                                                                      ece: ExecutionContextExecutor,
                                                                      materializer: ActorMaterializer)
  extends SprayJsonSupport with DefaultJsonProtocol with StatusCheckedSubsystem {

  override val statusUri = uri"$scheme://$interface:$port/engine/v1/status"

  def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val cromwellRequest = httpRequest
      .copy(uri = httpRequest.uri.withAuthority(interface, port).withScheme(scheme))
      .withHeaders(httpRequest.headers.filterNot(header => header.name == TimeoutAccessHeader))
    Http().singleRequest(cromwellRequest)
  } recoverWith {
    case e => Future.failed(CromwellConnectionFailure(e))
  }

  def extractWorkflowIds(response: HttpResponse, isBatch: Boolean): Future[List[String]] = {

    if (response.status.isSuccess) {
      if (isBatch) {
        Unmarshal(response).to[List[CromwellStatus]].map(_.map(_.id))
      } else {
        Unmarshal(response).to[CromwellStatus].map(s => List(s.id))
      }
    } else Future.successful(List.empty)
  }
}

object CromwellClient {
  // Header that Akka HTTP adds to every request on receive.
  // We get an warning in logs if we don't strip it out before sending the request to cromwell
  // HTTP header ‘Timeout-Access: <function1>’ is not allowed in requests
  // See: https://github.com/akka/akka-http/issues/64
  val TimeoutAccessHeader = "Timeout-Access"

  private case class CromwellConnectionFailure(f: Throwable) extends Exception(s"Unable to connect to Cromwell (${f.getMessage})", f)
}
