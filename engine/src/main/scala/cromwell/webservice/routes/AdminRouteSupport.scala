package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowAndState, WorkflowsBySubmissionId}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{ListSubmissions, ListSubmissionsResponse, ListSubmissionsResponseFailure, ListSubmissionsResponseSuccess}
import cromwell.services.admin.AdminServiceMessages._
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
  implicit val pauseSubmissionEncoder: Encoder[PauseSubmissionSuccess.type] = deriveEncoder[PauseSubmissionSuccess.type]

  val adminRoutes = concat(
    path("admin" / Segment / "listSubmissions") { _ =>
      onComplete(workflowStoreActor.ask(ListSubmissions).mapTo[ListSubmissionsResponse]) {
        case Success(response: ListSubmissionsResponseSuccess) =>
          complete(response.asJson)
        case Success(response: ListSubmissionsResponseFailure) =>
          new Exception(response.reason).failRequest(StatusCodes.BadRequest)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    },
    path("admin" / Segment / "pauseSubmission") { _ =>
      onComplete(serviceRegistryActor.ask(PauseSubmissionRequest).mapTo[PauseSubmissionResult]) {
        case Success(PauseSubmissionSuccess) =>
          complete(PauseSubmissionSuccess.asJson)
        case Success(response: PauseSubmissionFailure) =>
          new Exception(response.reason).failRequest(StatusCodes.BadRequest)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    }
  )
}
