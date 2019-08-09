package cromwell.webservice.routes

import java.util.UUID

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.{Multipart, StatusCodes}
import akka.http.scaladsl.server.Directives.{post, _}
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowAndState, WorkflowsBySubmissionId}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.services.metadata.MetadataService.{FetchWorkflowStatusFromMetadata, FetchWorkflowStatusInSummary, PutMetadataAction}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse, MetadataBuilderActorResponse}
import cromwell.webservice.routes.CromwellApiService.{InvalidWorkflowException, UnrecognizedWorkflowException}
import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
import io.circe.Encoder
import io.circe.generic.semiauto.deriveEncoder
import io.circe.syntax._
import spray.json.{JsObject, JsString}

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

trait AdminRouteSupport extends WebServiceUtils with MetadataRouteSupport {

  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef
  val workflowStoreActor: ActorRef

  implicit val workflowAndStateEncoder: Encoder[WorkflowAndState] = deriveEncoder[WorkflowAndState]
  implicit val listSubmissionsDataEncoder: Encoder[WorkflowsBySubmissionId] = deriveEncoder[WorkflowsBySubmissionId]
  implicit val listSubmissionsEncoder: Encoder[ListSubmissionsResponseSuccess] = deriveEncoder[ListSubmissionsResponseSuccess]

  implicit val pauseResponseEncoder: Encoder[PauseResponseSuccess] = deriveEncoder[PauseResponseSuccess]

  implicit val workflowStatusEncoder: Encoder[WorkflowStoreWorkflowStatus] = deriveEncoder[WorkflowStoreWorkflowStatus]
  implicit val fetchStatusResponseEncoder: Encoder[FetchWorkflowStatusResponseSuccess] = deriveEncoder[FetchWorkflowStatusResponseSuccess]

  val adminRoutes = concat(
    path("api"/ "admin" / Segment / "listSubmissions") { _ =>
      onComplete(workflowStoreActor.ask(ListSubmissions).mapTo[ListSubmissionsResponse]) {
        case Success(response: ListSubmissionsResponseSuccess) =>
          complete(response.asJson)
        case Success(response: ListSubmissionsResponseFailure) =>
          new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    },
    path("api"/ "admin" / Segment / "pauseSubmission") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          val pauseProcess = materializeFormData(formData) flatMap { args =>
            val submissionId = UUID.fromString(args("submissionId").utf8String)
            workflowStoreActor.ask(PauseSubmission(submissionId)).mapTo[PauseSubmissionResponse]
          }
          onComplete(pauseProcess) {
            case Success(resp: PauseResponseSuccess) =>
              complete(resp.asJson)
            case Success(response: PauseResponseFailure) =>
              new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
            case Failure(e) =>
              e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    },
    path("api"/ "admin" / Segment / "pauseAll") { _ =>
      post {
        val pauseProcess = workflowStoreActor.ask(PauseAll).mapTo[PauseSubmissionResponse]

        onComplete(pauseProcess) {
          case Success(resp: PauseResponseSuccess) =>
            complete(resp.asJson)
          case Success(response: PauseResponseFailure) =>
            new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
          case Failure(e) =>
            e.failRequest(StatusCodes.InternalServerError)
        }
      }
    },
    path("api"/ "admin" / Segment / "releaseHoldAcrossSubmission") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          val pauseProcess = materializeFormData(formData) flatMap { args =>
            val submissionId = UUID.fromString(args("submissionId").utf8String)
            val maxReleases = args.get("maxReleases").map(bytes => bytes.utf8String.toLong)
            workflowStoreActor.ask(WorkflowStoreActor.ReleaseHoldOnSubmission(submissionId, maxReleases)).mapTo[PauseSubmissionResponse]
          }
          onComplete(pauseProcess) {
            case Success(resp: PauseResponseSuccess) =>
              complete(resp.asJson)
            case Success(response: PauseResponseFailure) =>
              new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
            case Failure(e) =>
              e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    },
    path("admin" / Segment / Segment / "workflowStatuses") { (_, possibleWorkflowId) =>
      get {
        val allResponses = for {
          workflowId <- CromwellApiService.validateWorkflowIdInMetadata(possibleWorkflowId, serviceRegistryActor)
          responseA <- workflowStoreActor.ask(FetchWorkflowStatus(workflowId)).mapTo[FetchWorkflowStatusResponse]
          responseB <- metadataBuilderRegulatorActor.ask(FetchWorkflowStatusInSummary(workflowId)).mapTo[MetadataBuilderActorResponse]
          responseC <- metadataBuilderRegulatorActor.ask(FetchWorkflowStatusFromMetadata(workflowId)).mapTo[MetadataBuilderActorResponse]
        } yield (responseA, responseB, responseC)

        onComplete(allResponses) {
          case Success((a, b, c)) => {
            (a, b, c) match {
              case (r1: FetchWorkflowStatusResponseSuccess, r2: BuiltMetadataResponse, r3: BuiltMetadataResponse) => {
                val workflowStoreJsObject: JsObject = convertWorkflowStoreStatusToJsObject(r1.workflowStoreStatus)
                complete(combineJsObjects(workflowStoreJsObject, r2.response, r3.response))
              }
              case (r1: FetchWorkflowStatusResponseFailure, r2: FailedMetadataResponse, r3: FailedMetadataResponse) =>
                new Exception(r1.reason.getMessage + r2.reason.getMessage + r3.reason.getMessage).errorRequest(StatusCodes.InternalServerError)
              case (_, _, _) => new Exception("Something went wrong!").errorRequest(StatusCodes.InternalServerError)
            }
          }
          case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
          case Failure(e: InvalidWorkflowException) => e.failRequest(StatusCodes.BadRequest)
          case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
        }
      }
    },
    path("api"/ "admin" / Segment / "insertTerminalStatusInMetadata") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          val futureResult = materializeFormData(formData) map { args =>
            val workflowId = WorkflowId(UUID.fromString(args("id").utf8String))
            val terminalWorkflowStatus = WorkflowState.withName(args("terminalWorkflowStatus").utf8String)
            val event = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue(terminalWorkflowStatus))
            serviceRegistryActor ! PutMetadataAction(event)
            ()
          }

          onComplete(futureResult) {
            case Success(_) => complete(StatusCodes.OK)
            case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    }
  )

  private def convertWorkflowStoreStatusToJsObject(statusOption: Option[WorkflowStoreWorkflowStatus]): JsObject = {
    statusOption match {
      case Some(s) => JsObject(Map(
        WorkflowMetadataKeys.Status -> JsString(s.state),
        WorkflowMetadataKeys.SubmissionTime -> JsString(s.submissionTime.toString)
      ))
      case None => JsObject.empty
    }
  }

  private def combineJsObjects(js1: JsObject, js2: JsObject, js3: JsObject): JsObject = {
    JsObject(Map(
      "workflowStoreStatus" -> js1,
      "summaryTableStatus" -> js2,
      "metadataTableStatus" -> js3
    ))
  }
}
