package cromwell.webservice.routes.wes

import java.net.URL
import akka.actor.{ActorRef, ActorSystem}
import akka.http.scaladsl.model._
import akka.util.Timeout

import scala.concurrent.{ExecutionContext, Future}
import akka.stream.ActorMaterializer
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{WorkflowQueryFailure, WorkflowQueryResult, WorkflowQuerySuccess}
import cromwell.webservice.routes.MetadataRouteSupport.metadataQueryRequest
import cromwell.webservice.routes.wes.WesState._

final class Wes2CromwellInterface(cromwellPath: URL)(implicit system: ActorSystem, mat: ActorMaterializer, ec: ExecutionContext) {

  val serviceRegistryActor: ActorRef
  implicit val timeout: Timeout

  def listRuns(pageSize: Option[Int], pageToken: Option[String], headers: List[HttpHeader]): Future[WesResponse] = {
    // FIXME: to handle - page_size, page_token
    // FIXME: How to handle next_page_token in response?

    val metadataResponse = metadataQueryRequest(Seq.empty[(String, String)], serviceRegistryActor)
    metadataResponse.map {
      x => WesResponseRunList(fromMetadataQueryResponse(x)).runs
    }
  }

//    val cromwellUrl = s"$cromwellPath/query"
//    forwardToCromwell(cromwellUrl, headers, HttpMethods.GET , (s: String) => WesResponseRunList(RunListResponse.fromJson(s).runs))

  def fromMetadataQueryResponse(response: MetadataService.MetadataQueryResponse): RunListResponse = {
    response match {
      case w: WorkflowQuerySuccess =>
        val runs = w.response.results.toList.map(x => WesRunStatus(x.id, fromABC(x.status)))
        RunListResponse(runs, "Not Yet Implemented")
      case w: WorkflowQueryFailure => ???
    }
}

  def listRuns(pageSize: Option[Int], pageToken: Option[String], headers: List[HttpHeader]): Future[WesResponse] = {
    // FIXME: to handle - page_size, page_token
    // FIXME: How to handle next_page_token in response?
    val cromwellUrl = s"$cromwellPath/query"
    forwardToCromwell(cromwellUrl, headers, HttpMethods.GET , (s: String) => WesResponseRunList(RunListResponse.fromJson(s).runs))
  }