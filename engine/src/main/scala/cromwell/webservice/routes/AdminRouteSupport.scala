package cromwell.webservice.routes

import java.util.UUID

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.model.{Multipart, StatusCodes}
import akka.http.scaladsl.server.Directives.{post, _}
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.WorkflowSubmitted
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowAndState, WorkflowsBySubmissionId}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
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
      }
  )
}
