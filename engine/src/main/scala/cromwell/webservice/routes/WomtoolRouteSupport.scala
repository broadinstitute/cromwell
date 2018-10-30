package cromwell.webservice.routes

import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._

import scala.concurrent.{ExecutionContext, Future}

trait WomtoolRouteSupport extends {

  implicit val ec: ExecutionContext

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        entity(as[Multipart.FormData]) { _: Multipart.FormData =>
          complete(Future(ToResponseMarshallable("A")))
        }
      }
    }

}
