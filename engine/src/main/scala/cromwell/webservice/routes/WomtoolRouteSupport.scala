package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.{WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeFailure, DescribeRequest, DescribeResult, DescribeSuccess}
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

trait WomtoolRouteSupport extends WebServiceUtils {

  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          onComplete(materializeFormData(formData)) {
            case Success(data) =>
              validateAndSubmitRequest(data)
            case Failure(e) =>
              e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    }

  private def validateAndSubmitRequest(data: MaterializedFormData): Route = {
    // TODO: move constants to WebServiceUtils, adopt in PartialWorkflowSources
    val workflowSource = data.get("workflowSource").map(_.utf8String)
    val workflowUrl = data.get("workflowUrl").map(_.utf8String)
    val workflowInputs = data.get("workflowInputs").map(_.utf8String)
    val workflowType = data.get("workflowType").map(_.utf8String)
    val workflowVersion = data.get("workflowTypeVersion").map(_.utf8String)

    val wsfc = WorkflowSourceFilesCollection(
      workflowSource,
      workflowUrl,
      workflowRoot = None,
      workflowType,
      workflowVersion,
      workflowInputs.getOrElse(""),
      workflowOptions = WorkflowOptions.empty,
      labelsJson = "",
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )

    onComplete(serviceRegistryActor.ask(DescribeRequest(wsfc)).mapTo[DescribeResult]) {
      case Success(response: DescribeSuccess) =>
        import cromwell.services.womtool.models.WorkflowDescription.workflowDescriptionEncoder
        import de.heikoseeberger.akkahttpcirce.ErrorAccumulatingCirceSupport._
        import io.circe.syntax._

        complete(response.description.asJson)
      case Success(response: DescribeFailure) =>
        new Exception(response.reason).failRequest(StatusCodes.BadRequest)
      case Failure(e) =>
        e.failRequest(StatusCodes.InternalServerError)
    }
  }

}
