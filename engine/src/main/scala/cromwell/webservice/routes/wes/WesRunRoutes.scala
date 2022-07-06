package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.{StatusCode, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.{Directive1, Route}
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.pattern.AskTimeoutException
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
// import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.WorkflowStoreSubmitActorResponse
import cromwell.server.CromwellShutdown
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, GetSingleWorkflowMetadataAction}
import cromwell.services.{FailedMetadataJsonResponse, SuccessfulMetadataJsonResponse}
import cromwell.webservice.routes.CromwellApiService
import cromwell.webservice.routes.MetadataRouteSupport.{metadataBuilderActorRequest, metadataQueryRequest}
import cromwell.webservice.routes.wes.WesResponseJsonSupport.{WesResponseErrorFormat, WesResponseFormat}
import cromwell.webservice.routes.wes.WesRunRoutes.{completeCromwellResponse, extractSubmission, runLog}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success}

trait WesRunRoutes extends CromwellApiService {

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
            post {
              extractSubmission() { submission =>
                submitRequest(submission.entity,
                  isSingleSubmission = true,
                  //successHandler = WesSuccessHandler,
                  //errorHandler = WesErrorHandler
                )
              }
            }
        },
        path(Segment) { workflowId =>
          get {
            // this is what it was like in code found in the project… it perhaps isn’t ideal but doesn’t seem to hurt, so leaving it like this for now.
            completeCromwellResponse(runLog(workflowId, (w: WorkflowId) => GetSingleWorkflowMetadataAction(w, None, None, expandSubWorkflows = false), serviceRegistryActor))
          }
      }
      )
    }
}

object WesRunRoutes {

  implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

//  def WesSuccessHandler: PartialFunction[WorkflowStoreSubmitActorResponse, Route] = {
//    case response => complete(WesRunId(response.workflowId.toString))
//  }

  def WesErrorHandler: PartialFunction[Throwable, Route] = {
    case e: IllegalStateException => respondWithWesError(e.getLocalizedMessage, StatusCodes.Forbidden)
    case e: WorkflowNotFoundException => respondWithWesError(e.getLocalizedMessage, StatusCodes.NotFound)
    case _: AskTimeoutException if CromwellShutdown.shutdownInProgress() => respondWithWesError("Cromwell service is shutting down", StatusCodes.InternalServerError)
    case e: Exception => respondWithWesError(e.getLocalizedMessage, StatusCodes.InternalServerError)
  }

  private def respondWithWesError(errorMsg: String, status: StatusCode): Route = {
    complete((status, WesErrorResponse(errorMsg, status.intValue)))
  }


  def extractSubmission(): Directive1[WesSubmission] = {
    formFields((
      "workflow_params".?,
      "workflow_type".?,
      "workflow_type_version".?,
      "tags".?,
      "workflow_engine_parameters".?,
      "workflow_url".?,
      "workflow_attachment".as[String].*
    )).as(WesSubmission)

    }

  def completeCromwellResponse(future: => Future[WesResponse]): Route = {
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

    metadataJsonResponse.map {
      case SuccessfulMetadataJsonResponse(_, responseJson) => WesResponseWorkflowMetadata(WesRunLog.fromJson(responseJson.toString()))
      case FailedMetadataJsonResponse(_, reason) => WesErrorResponse(reason.getMessage, StatusCodes.InternalServerError.intValue)
    }
  }
}
