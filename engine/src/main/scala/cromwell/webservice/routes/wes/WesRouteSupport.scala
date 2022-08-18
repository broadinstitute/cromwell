package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.model.{StatusCode, StatusCodes}
import akka.http.scaladsl.server.Directives.{path, _}
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.http.scaladsl.server.{Directive1, Route}
import akka.pattern.{AskTimeoutException, ask}
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.core.abort.SuccessfulAbortResponse
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.server.CromwellShutdown
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, GetSingleWorkflowMetadataAction, GetStatus, MetadataServiceResponse, StatusLookupFailed}
import cromwell.services.{FailedMetadataJsonResponse, SuccessfulMetadataJsonResponse}
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import cromwell.webservice.routes.CromwellApiService.{UnrecognizedWorkflowException, validateWorkflowIdInMetadata}
import cromwell.webservice.routes.MetadataRouteSupport.{metadataBuilderActorRequest, metadataQueryRequest}
import cromwell.webservice.routes.wes.WesResponseJsonSupport._
import cromwell.webservice.routes.wes.WesRouteSupport._
import cromwell.webservice.routes.{CromwellApiService, WesCromwellRouteSupport}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}


trait WesRouteSupport extends HttpInstrumentation with WesCromwellRouteSupport {
  val serviceRegistryActor: ActorRef
  val workflowManagerActor: ActorRef
  val workflowStoreActor: ActorRef

  implicit val ec: ExecutionContext
  implicit val timeout: Timeout

  /*
    Defines routes intended to sit alongside the primary Cromwell REST endpoints. For instance, we'll now have:
    GET /api/workflows/v1/ID/status
    GET /ga4gh/wes/v1/runs/ID

    POST /api/workflows/v1/ID/abort
    POST /ga4gh/wes/v1/runs/ID/cancel

    These routes are not currently going through the MetadataBuilderRegulator/MetadataBuilder. This is for two reasons:
      - It'd require a fairly substantial refactor of the MetadataBuilderActor to be more general
      - It's expected that for now the usage of these endpoints will not be extensive, so the protections of the regulator
         should not be necessary
    */
  val wesRoutes: Route =
    instrumentRequest {
      concat(
        pathPrefix("ga4gh" / "wes" / "v1") {
          concat(
            path("service-info") {
              get {
                complete(ServiceInfo.toWesResponse(workflowStoreActor))
              }
            },
            path("runs") {
              get {
                parameters(("page_size".as[Int].?, "page_token".?)) { (pageSize, pageToken) =>
                  completeCromwellResponse(listRuns(pageSize, pageToken, serviceRegistryActor))
                }
              } ~
                post {
                  extractSubmission() { submission =>
                    submitRequest(submission.entity,
                      isSingleSubmission = true)
                  }
                }
            },
            path("runs" / Segment) { workflowId =>
              get {
                // this is what it was like in code found in the project… it perhaps isn’t ideal but doesn’t seem to hurt, so leaving it like this for now.
                completeCromwellResponse(runLog(workflowId, (w: WorkflowId) => GetSingleWorkflowMetadataAction(w, None, None, expandSubWorkflows = false), serviceRegistryActor))
              }
            },
              path("runs" / Segment / "status") { possibleWorkflowId =>
                val response = validateWorkflowIdInMetadata(possibleWorkflowId, serviceRegistryActor).flatMap(w => serviceRegistryActor.ask(GetStatus(w)).mapTo[MetadataServiceResponse])
                // WES can also return a 401 or a 403 but that requires user auth knowledge which Cromwell doesn't currently have
                onComplete(response) {
                  case Success(SuccessfulMetadataJsonResponse(_, jsObject)) =>
                    val wesState = WesState.fromCromwellStatusJson(jsObject)
                    complete(WesRunStatus(possibleWorkflowId, wesState))
                  case Success(r: StatusLookupFailed) => r.reason.errorRequest(StatusCodes.InternalServerError)
                  case Success(m: MetadataServiceResponse) =>
                    // This should never happen, but ....
                    val error = new IllegalStateException("Unexpected response from Metadata service: " + m)
                    error.errorRequest(StatusCodes.InternalServerError)
                  case Failure(_: UnrecognizedWorkflowException) => complete(NotFoundError)
                  case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
                }
              },
            path("runs" / Segment / "cancel") { possibleWorkflowId =>
              post {
                CromwellApiService.abortWorkflow(possibleWorkflowId,
                  workflowStoreActor,
                  workflowManagerActor,
                  successHandler = WesAbortSuccessHandler,
                  errorHandler = WesAbortErrorHandler)
              }
            }
          )
        }
      )
    }
}

object WesRouteSupport {
  import WesResponseJsonSupport._

  implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

  val NotFoundError = WesErrorResponse("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)

  def WesAbortSuccessHandler: PartialFunction[SuccessfulAbortResponse, Route] = {
    case response => complete(WesRunId(response.workflowId.toString))
  }

  def WesAbortErrorHandler: PartialFunction[Throwable, Route] = {
    // There are also some auth situations which should be handled, but at the moment Cromwell doesn't allow for those
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
