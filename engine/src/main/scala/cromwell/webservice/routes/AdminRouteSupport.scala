package cromwell.webservice.routes

import java.util.UUID

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.model.{Multipart, StatusCodes}
import akka.http.scaladsl.server.Directives.{post, _}
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState, WorkflowSubmitted}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowAndState, WorkflowsBySubmissionId}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
import io.circe.Encoder
import io.circe.generic.semiauto.deriveEncoder
import io.circe.syntax._

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

trait AdminRouteSupport extends WebServiceUtils {

  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef
  val workflowStoreActor: ActorRef

  implicit val workflowAndStateEncoder: Encoder[WorkflowAndState] = deriveEncoder[WorkflowAndState]
  implicit val listSubmissionsDataEncoder: Encoder[WorkflowsBySubmissionId] = deriveEncoder[WorkflowsBySubmissionId]
  implicit val listSubmissionsEncoder: Encoder[ListSubmissionsResponseSuccess] = deriveEncoder[ListSubmissionsResponseSuccess]

  implicit val pauseResponseEncoder: Encoder[PauseSubmissionResponseSuccess] = deriveEncoder[PauseSubmissionResponseSuccess]

  val adminRoutes = concat(
    path("admin" / Segment / "listSubmissions") { _ =>
      onComplete(workflowStoreActor.ask(ListSubmissions).mapTo[ListSubmissionsResponse]) {
        case Success(response: ListSubmissionsResponseSuccess) =>
          complete(response.asJson)
        case Success(response: ListSubmissionsResponseFailure) =>
          new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    },
    path("admin" / Segment / "pauseSubmission") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          onComplete(materializeFormData(formData) flatMap { args =>
            val submissionId = UUID.fromString(args("submissionId").utf8String)

            workflowStoreActor.ask(PauseSubmission(submissionId, WorkflowSubmitted)).mapTo[PauseSubmissionResponse] }) {
              case Success(resp: PauseSubmissionResponseSuccess) =>
                complete(resp.asJson)
              case Success(response: PauseSubmissionResponseFailure) =>
                new Exception(response.reason).failRequest(StatusCodes.InternalServerError)
              case Failure(e) =>
                e.failRequest(StatusCodes.InternalServerError)
            }
          }
        }
      },
    path("admin" / Segment / "insertTerminalStatusInMetadata") { _ =>
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
}
