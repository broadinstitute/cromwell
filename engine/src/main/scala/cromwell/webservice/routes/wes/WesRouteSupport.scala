package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.{AskTimeoutException, ask}
import akka.util.Timeout
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.services.metadata.MetadataService.{GetSingleWorkflowMetadataAction, GetStatus, MetadataServiceResponse, StatusLookupFailed, StatusLookupResponse}
import cromwell.webservice.routes.CromwellApiService.{UnrecognizedWorkflowException, validateWorkflowId}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}
import WesResponseJsonSupport._
import akka.http.scaladsl.model.{StatusCode, StatusCodes}
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import WesRouteSupport._
import cats.data.NonEmptyList
import cromwell.core.WorkflowId
import cromwell.core.abort.SuccessfulAbortResponse
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.server.CromwellShutdown
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse}
import cromwell.webservice.routes.{CromwellApiService, MetadataRouteSupport}

trait WesRouteSupport extends HttpInstrumentation {
  val serviceRegistryActor: ActorRef
  val workflowManagerActor: ActorRef
  val workflowStoreActor: ActorRef
  val metadataBuilderRegulatorActor: ActorRef

  implicit val ec: ExecutionContext
  implicit val timeout: Timeout

  /*
    Defines routes intended to sit alongside the primary Cromwell REST endpoints. For instance, we'll now have:
    GET /api/workflows/v1/ID/status
    GET /ga4gh/wes/v1/runs/ID

    POST /api/workflows/v1/ID/abort
    POST /ga4gh/wes/v1/runs/ID/cancel

    Most of these routes are not currently going through the MetadataBuilderRegulator/MetadataBuilder. This is for two reasons:
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
                path(Segment) { workflowId =>
                  get {
                    /*
                  This is the one exception to the "not going through the MetadataBuilderRegulator/MetadataBuilder"
                  statement above. This endpoint is basically just reshaping the metadata altogether so it still
                  needs to be built out of the DB, but MetadataBuilder starts forming its json all the way at the
                  start of the process. A longer term project would be to refactor the building process to be more flexible
                  but for now we'll grab the built metadata, parse it back into scala objects of the right shape and
                  then push it back out the other side
                 */
                    val metadata = MetadataRouteSupport.metadataBuilderActorRequest(workflowId,
                      (w: WorkflowId) => GetSingleWorkflowMetadataAction(WorkflowId.fromString(workflowId),
                        MetadataIncludeKeys,
                        None,
                        true), serviceRegistryActor, metadataBuilderRegulatorActor)
                    onComplete(metadata) {
                      case Success(r: BuiltMetadataResponse) => complete(WesResponseRunLog.fromJson(r.response))
                      case Success(f: FailedMetadataResponse) => respondWithWesError(f.reason.getMessage, StatusCodes.InternalServerError.intValue)
                      case Failure(_: UnrecognizedWorkflowException) => respondWithNotFound
                      case Failure(f) => respondWithWesError(f.getMessage, StatusCodes.InternalServerError.intValue)
                    }
                  }
                },
                path(Segment / "status") { possibleWorkflowId =>
                  val response = validateWorkflowId(possibleWorkflowId, serviceRegistryActor).flatMap(w => serviceRegistryActor.ask(GetStatus(w)).mapTo[MetadataServiceResponse])
                  // WES can also return a 401 or a 403 but that requires user auth knowledge which Cromwell doesn't currently have
                  onComplete(response) {
                    case Success(s: StatusLookupResponse) =>
                      val wesState = WesState.fromCromwellStatus(s.status)
                      complete(WesRunStatus(s.workflowId.toString, wesState))
                    case Success(r: StatusLookupFailed) => respondWithWesError(r.reason.getMessage, StatusCodes.InternalServerError)
                    case Success(m: MetadataServiceResponse) =>
                      // This should never happen, but ....
                      respondWithWesError("Unexpected response from Metadata service: " + m, StatusCodes.InternalServerError)
                    case Failure(_: UnrecognizedWorkflowException) => respondWithNotFound
                    case Failure(e) => respondWithWesError(e.getMessage, StatusCodes.InternalServerError.intValue)
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
  def notFoundError = respondWithWesError("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)

  val MetadataIncludeKeys = Option(NonEmptyList.of("shardIndex", "commandLine", "returnCode", "start", "end", "stdout",
                                            "stderr", "workflow", "options", "inputs", "labels", "workflowName", "id",
                                            "status", "submittedFiles", "outputs", "actualWorkflowLanguage", "ActualWorkflowLanguageVersion",
                                             "workflowUrl"))

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

  def respondWithNotFound = respondWithWesError("The requested workflow run wasn't found", StatusCodes.NotFound.intValue)

  def respondWithWesError(errorMsg: String, status: StatusCode): Route = {
    complete((status, WesErrorResponse(errorMsg, status.intValue)))
  }
}
