package cromwell.webservice.routes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.webservice.WebServiceUtils
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import spray.json.{JsBoolean, JsObject}
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

trait WomtoolRouteSupport extends WebServiceUtils {

  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          onComplete(materializeFormData(formData)) {
            case Success(data) =>

              // TODO: move constants to FormDataSupport, adopt in PartialWorkflowSources
              val workflowSource = data.get("workflowSource").map(_.utf8String)
              val workflowUrl = data.get("workflowUrl").map(_.utf8String)

              (workflowSource, workflowUrl) match {
                case (Some(_), Some(_)) =>
                  new IllegalArgumentException("Must submit exactly one of workflow source, workflow URL; received both").failRequest(StatusCodes.BadRequest)
                case (None, Some(_)) =>
                  new Exception("URL submissions not yet supported").failRequest(StatusCodes.NotImplemented)
                case (Some(source), None) =>
                  complete {
                    // After much discussion, it was resolved that returning Future(something) is the right way to do
                    // things in Akka HTTP, when no legacy or migration constraints funnel one into using PerRequest
                    doExpensiveWDLValidation(source) map { result: Boolean =>
                      JsObject(Map("valid" -> JsBoolean(result)))
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

  private def doExpensiveWDLValidation(workflow: String): Future[Boolean] = Future(math.random < 0.5)

}
