package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.util.ImportResolver.HttpResolver
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeRequest, DescribeResponse}
import cromwell.services.womtool.WomtoolServiceMessages.JsonSupport.describeResponseFormat
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import wom.core.WorkflowSource

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
      workflowOptionsJson = "",
      labelsJson = "",
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )

    // The HTTP resolver is used to pull down workflows submitted by URL
    // TODO: we are potentially tying up the web thread here making a request to the workflow URL?
    LanguageFactoryUtil.findWorkflowSource(workflowSource, workflowUrl, List(HttpResolver(None, Map.empty))) match {
      case Right(sourceAndResolvers: (WorkflowSource, List[ImportResolver.ImportResolver])) =>
        complete {
          serviceRegistryActor.ask(DescribeRequest(sourceAndResolvers._1, wsfc)) map {
            case result: DescribeResponse =>
              result
          }
        }
      case Left(errors) =>
        new Exception(errors.toList.mkString(", ")).failRequest(StatusCodes.BadRequest)
    }
  }

}
