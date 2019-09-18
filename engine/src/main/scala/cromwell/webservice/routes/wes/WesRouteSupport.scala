package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.{AskTimeoutException, ask}
import akka.util.Timeout
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.services.metadata.MetadataService.{GetStatus, MetadataServiceResponse, StatusLookupFailed}
import cromwell.webservice.routes.CromwellApiService.{UnrecognizedWorkflowException, validateWorkflowIdInMetadata}
import cromwell.webservice.WebServiceUtils.EnhancedThrowable

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}
import WesResponseJsonSupport._
import akka.http.scaladsl.model.{StatusCode, StatusCodes}
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import WesRouteSupport._
import cromwell.core.abort.SuccessfulAbortResponse
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.server.CromwellShutdown
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.BuiltMetadataResponse
import cromwell.webservice.routes.CromwellApiService

trait WesRouteSupport extends HttpInstrumentation {
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
              complete(ServiceInfo.toWesResponse(workflowStoreActor))
            },
            pathPrefix("runs") {
              concat(
                path(Segment / "status") { possibleWorkflowId =>
                  val response = validateWorkflowIdInMetadata(possibleWorkflowId, serviceRegistryActor).flatMap(w => serviceRegistryActor.ask(GetStatus(w)).mapTo[MetadataServiceResponse])
                  // WES can also return a 401 or a 403 but that requires user auth knowledge which Cromwell doesn't currently have
                  onComplete(response) {
                    case Success(BuiltMetadataResponse(_, jsObject)) =>
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
                path(Segment / "cancel") { possibleWorkflowId =>
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
      )
    }
}

object WesRouteSupport {
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
}
