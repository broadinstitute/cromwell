package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.model.{Multipart, StatusCode, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.http.scaladsl.server.{Directive1, Route}
import akka.pattern.{AskTimeoutException, ask}
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.abort.SuccessfulAbortResponse
import cromwell.core.{WorkflowId, WorkflowOnHold, WorkflowState, WorkflowSubmitted}
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreSubmitActor}
import cromwell.server.CromwellShutdown
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, GetSingleWorkflowMetadataAction, GetStatus, MetadataServiceResponse, StatusLookupFailed}
import cromwell.services.{FailedMetadataJsonResponse, SuccessfulMetadataJsonResponse}
import cromwell.webservice.PartialWorkflowSources
import cromwell.webservice.WebServiceUtils.{EnhancedThrowable, completeResponse, materializeFormData}
import cromwell.webservice.routes.CromwellApiService
import cromwell.webservice.routes.CromwellApiService.{UnrecognizedWorkflowException, validateWorkflowIdInMetadata}
import cromwell.webservice.routes.MetadataRouteSupport.{metadataBuilderActorRequest, metadataQueryRequest}
import cromwell.webservice.routes.wes.WesResponseJsonSupport._
import cromwell.webservice.routes.wes.WesRouteSupport.{respondWithWesError, _}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future, TimeoutException}
import scala.util.{Failure, Success}



trait WesRouteSupport extends HttpInstrumentation {

  val serviceRegistryActor: ActorRef
  val workflowManagerActor: ActorRef
  val workflowStoreActor: ActorRef

  implicit val ec: ExecutionContext
  implicit val timeout: Timeout
  implicit val materializer: ActorMaterializer

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
                    wesSubmitRequest(submission.entity,
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

  def toWesResponse(workflowId: WorkflowId, workflowState: WorkflowState):  WesRunStatus = {
    WesRunStatus(workflowId.toString, WesState.fromCromwellStatus(workflowState))
  }

  def toWesResponseId(workflowId: WorkflowId): WesRunId ={
    WesRunId(workflowId.toString)
  }

  def wesSubmitRequest(formData: Multipart.FormData, isSingleSubmission: Boolean): Route = {
    def getWorkflowState(workflowOnHold: Boolean): WorkflowState = {
      if (workflowOnHold)
        WorkflowOnHold
      else WorkflowSubmitted
    }

    def sendToWorkflowStore(command: WorkflowStoreActor.WorkflowStoreActorSubmitCommand, warnings: Seq[String], workflowState: WorkflowState): Route = {
      // NOTE: Do not blindly copy the akka-http -to- ask-actor pattern below without knowing the pros and cons.
      onComplete(workflowStoreActor.ask(command).mapTo[WorkflowStoreSubmitActor.WorkflowStoreSubmitActorResponse]) {
        case Success(w) =>
          w match {
            case WorkflowStoreSubmitActor.WorkflowSubmittedToStore(workflowId, _) =>
              completeResponse(StatusCodes.Created, toWesResponseId(workflowId), warnings)
            case WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore(workflowIds, _) =>
              completeResponse(StatusCodes.Created, workflowIds.toList.map(toWesResponse(_, workflowState)), warnings)
            case WorkflowStoreSubmitActor.WorkflowSubmitFailed(throwable) =>
              respondWithWesError(throwable.getLocalizedMessage, StatusCodes.BadRequest)
          }
        case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => respondWithWesError("Cromwell service is shutting down", StatusCodes.InternalServerError)
        case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
        case Failure(e) => e.failRequest(StatusCodes.InternalServerError, warnings)
      }
    }

    onComplete(materializeFormData(formData)) {
      case Success(data) =>
        PartialWorkflowSources.fromSubmitRoute(data, allowNoInputs = isSingleSubmission) match {
          case Success(workflowSourceFiles) if isSingleSubmission && workflowSourceFiles.size == 1 =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            sendToWorkflowStore(WorkflowStoreActor.SubmitWorkflow(workflowSourceFiles.head), warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          // Catches the case where someone has gone through the single submission endpoint w/ more than one workflow
          case Success(workflowSourceFiles) if isSingleSubmission =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            val e = new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint.")
            e.failRequest(StatusCodes.BadRequest, warnings)
          case Success(workflowSourceFiles) =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            sendToWorkflowStore(
              WorkflowStoreActor.BatchSubmitWorkflows(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList)),
              warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          case Failure(t) => t.failRequest(StatusCodes.BadRequest)
        }
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => respondWithWesError(e.getLocalizedMessage, StatusCodes.InternalServerError)
    }
  }
}



object WesRouteSupport {
  import WesResponseJsonSupport._

  implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration
  import scala.concurrent.ExecutionContext.Implicits.global

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
      case SuccessfulMetadataJsonResponse(_, responseJson) => WesRunLog.fromJson(responseJson.toString())
      case FailedMetadataJsonResponse(_, reason) => WesErrorResponse(reason.getMessage, StatusCodes.InternalServerError.intValue)
    }
  }
}