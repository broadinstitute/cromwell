package cromwell.webservice.routes.wes

import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.model.HttpHeader
import akka.stream.ActorMaterializer
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import cromwell.services.metadata.MetadataService
import cromwell.webservice.routes.MetadataRouteSupport.metadataQueryRequest

import java.net.URL
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

final class Wes2CromwellInterface()(implicit system: ActorSystem, mat: ActorMaterializer, ec: ExecutionContext) {


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
