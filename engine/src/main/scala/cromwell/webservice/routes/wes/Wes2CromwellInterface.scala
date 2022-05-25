package cromwell.webservice.routes.wes

import akka.actor.{ActorRef} // ActorSystem
import akka.http.scaladsl.model.{HttpHeader, StatusCodes}
// import akka.stream.ActorMaterializer
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import cromwell.services.metadata.MetadataService
import cromwell.webservice.routes.MetadataRouteSupport.metadataQueryRequest

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

final class Wes2CromwellInterface()(implicit ec: ExecutionContext) {


  implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

  def listRuns(pageSize: Option[Int], pageToken: Option[String], headers: List[HttpHeader], serviceRegistryActor: ActorRef): Future[WesResponse] = {
    // FIXME: to handle - page_size, page_token
    // FIXME: How to handle next_page_token in response?
    val metadataResponse = metadataQueryRequest(Seq.empty[(String, String)], serviceRegistryActor)

    metadataResponse.map {
      x: MetadataService.MetadataQueryResponse => WesResponseRunList(RunListResponse.fromMetadataQueryResponse(x).runs)

    }
  }
}

object Wes2CromwellInterface {
  val BadRequestError = WesErrorResponse("The request is malformed", StatusCodes.BadRequest.intValue)
  val InternalServerError = WesErrorResponse("Cromwell server error", StatusCodes.InternalServerError.intValue)
  val UnauthorizedError = WesErrorResponse("The request is unauthorized", StatusCodes.Unauthorized.intValue)
  val ForbiddenError = WesErrorResponse("The requester is not authorized to perform this action", StatusCodes.Forbidden.intValue)
  val NotFoundError = WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)
}
