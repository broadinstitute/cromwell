package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.pattern.ask
import akka.util.Timeout
import cromwell.services.admin.AdminServiceMessages.{AdminResult, ListSubmissionsRequest, PauseSubmissionSuccess}
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

trait AdminRouteSupport extends WebServiceUtils {

  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef

  val adminRoutes = concat(
    path("admin" / Segment / "listSubmission") { _ =>
      onComplete(serviceRegistryActor.ask(ListSubmissionsRequest).mapTo[AdminResult]) {
        case Success(PauseSubmissionSuccess) =>
          import cromwell.services.womtool.models.WorkflowDescription.workflowDescriptionEncoder
          import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
          import io.circe.syntax._

          complete(PauseSubmissionSuccess.toString)
//        case Success(response: DescribeFailure) => ??? fill it in
//          new Exception(response.reason).failRequest(StatusCodes.BadRequest)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    },
    path("admin" / Segment / "pauseSubmission") { _ =>
      onComplete(serviceRegistryActor.ask(ListSubmissionsRequest).mapTo[AdminResult]) {
        case Success(PauseSubmissionSuccess) =>
          import cromwell.services.womtool.models.WorkflowDescription.workflowDescriptionEncoder
          import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
          import io.circe.syntax._

          complete(PauseSubmissionSuccess.toString)
        //        case Success(response: DescribeFailure) => ??? fill it in
        //          new Exception(response.reason).failRequest(StatusCodes.BadRequest)
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError)
      }
    }
  )
}
