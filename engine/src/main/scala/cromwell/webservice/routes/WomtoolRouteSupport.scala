package cromwell.webservice.routes

import akka.actor.ActorRef
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeRequest, DescribeResponse}
import cromwell.services.womtool.WomtoolServiceMessages.JsonSupport.describeResponseFormat
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

trait WomtoolRouteSupport extends WebServiceUtils {

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

              // TODO: would be nice to reuse findWorkflowSource
              (workflowSource, workflowUrl) match {
                case (Some(_), Some(_)) =>
                  new IllegalArgumentException("Must submit exactly one of workflow source, workflow URL; received both").failRequest(StatusCodes.BadRequest)
                case (None, Some(_)) =>
                  new Exception("URL submissions not yet supported").failRequest(StatusCodes.NotImplemented)
                case (Some(source), None) =>
                  complete {
                    // !!! The following is no longer true after service registry integration !!!
                    // !!! Perhaps another discussions is required? !!!

                    // After much discussion, it was resolved that returning Future(something) is the right way to do
                    // things in Akka HTTP, when no legacy or migration constraints funnel one into using PerRequest

                    serviceRegistryActor.ask(DescribeRequest(source, wsfc)) map {
                      case result: DescribeResponse =>
                        result
                    }
                  }
                case (None, None) =>
                  new IllegalArgumentException("Must submit exactly one of workflow source, workflow URL; received neither").failRequest(StatusCodes.BadRequest)
              }
            case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
          }
        }
      }
    }

}
