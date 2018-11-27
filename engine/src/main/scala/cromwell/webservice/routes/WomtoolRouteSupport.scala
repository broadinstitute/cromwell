package cromwell.webservice.routes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import spray.json.{JsObject, JsString}

import scala.concurrent.{ExecutionContext, Future}

trait WomtoolRouteSupport extends {

  implicit val ec: ExecutionContext

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { _: Multipart.FormData =>

          // After much discussion, it was resolved that returning Future(something) is the right way to do
          // things in Akka HTTP, when no legacy or migration constraints funnel one into using PerRequest
          complete(Future(JsObject(Map("a" -> JsString("b")))))
        }
      }
    }

}
