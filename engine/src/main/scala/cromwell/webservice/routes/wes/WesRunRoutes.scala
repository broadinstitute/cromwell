package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.services.SuccessfulMetadataJsonResponse
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, GetSingleWorkflowMetadataAction}
import cromwell.webservice.routes.MetadataRouteSupport.{metadataBuilderActorRequest, metadataQueryRequest}

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Failure, Success}

trait WesRunRoutes {

  val serviceRegistryActor: ActorRef

  lazy val runRoutes: Route =
    pathPrefix("ga4gh" / "wes" / "v1") {
      concat(
        pathPrefix("runs") {
          get {
            parameters(("page_size".as[Int].?, "page_token".?)) { (pageSize, pageToken) =>
              WesRunRoutes.completeCromwellResponse(WesRunRoutes.listRuns(pageSize, pageToken, serviceRegistryActor))
            }
          }
          concat(
            path(Segment) { workflowId =>
              get {
                get {
                  completeCromwellResponse(runLog(workflowId))
              }
            }
          )
        }
      )
    }
}

object WesRunRoutes {

import akka.util.Timeout
import scala.concurrent.duration.FiniteDuration
import net.ceedubs.ficus.Ficus._

implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
implicit lazy val timeout: Timeout = duration

def completeCromwellResponse(future: => Future[WesResponse]): Route = {

  import WesResponseJsonSupport.WesResponseErrorFormat
  import cromwell.webservice.routes.wes.WesResponseJsonSupport.WesResponseFormat
  import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._

  onComplete(future) {
    case Success(response: WesResponse) => complete(response)
    case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
  }
}

def listRuns(pageSize: Option[Int], pageToken: Option[String], serviceRegistryActor: ActorRef): Future[WesResponse] = {
  // FIXME: to handle - page_size, page_token
  // FIXME: How to handle next_page_token in response?
  metadataQueryRequest(Seq.empty[(String, String)], serviceRegistryActor).map(RunListResponse.fromMetadataQueryResponse)
  }

def runLog(workflowId: String, request: WorkflowId => BuildMetadataJsonAction, serviceRegistryActor: ActorRef): Future[WesResponse] = {

   val metadataJsonResponse = metadataBuilderActorRequest(workflowId, request, serviceRegistryActor)



  def fromCromwellMetadata(response: Future[MetadataJsonResponse], workflowId: String): WesResponse = {

    response match {
      case w: SuccessfulMetadataJsonResponse =>
        val runs = w.

      case e: FailedMetadataJsonResponse =>
    }
  }
}



  WesResponseWorkflowMetadata(WesRunLog.fromCromwellMetadata())
}

}
