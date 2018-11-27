package cromwell.webservice.routes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.webservice.FormDataSupport
import spray.json.{JsObject, JsString}

import scala.concurrent.ExecutionContext

trait WomtoolRouteSupport extends FormDataSupport {

  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { formData: Multipart.FormData =>
          complete {
            // After much discussion, it was resolved that returning Future(something) is the right way to do
            // things in Akka HTTP, when no legacy or migration constraints funnel one into using PerRequest
            materializeFormData(formData) map { formData: MaterializedFormData =>
              JsObject(Map("form size" -> JsString(formData.size.toString)))
            }
          }
        }
      }
    }

}
