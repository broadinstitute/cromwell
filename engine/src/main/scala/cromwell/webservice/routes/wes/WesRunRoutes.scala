package cromwell.webservice.routes.wes

import akka.actor.ActorRef
import akka.http.scaladsl.model.{HttpHeader, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.util.Timeout
import WesRunRoutes._
import akka.stream.ActorMaterializer

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Failure, Success}

trait WesRunRoutes {
  implicit def materializer: ActorMaterializer

  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef

  lazy val wes2CromwellInterface = new Wes2CromwellInterface()

  lazy val runRoutes: Route =
    optionalHeaderValue(extractAuthorizationHeader) { authHeader =>
      val cromwellRequestHeaders = authHeader.toList
      pathPrefix("ga4gh" / "wes" / "v1") {
        concat(
          pathPrefix("runs") {
            concat(
              pathEnd {
                concat(
                  get {
                    parameters(("page_size".as[Int].?, "page_token".?)) { (pageSize, pageToken) =>
                      completeCromwellResponse(wes2CromwellInterface.listRuns(pageSize, pageToken, cromwellRequestHeaders, serviceRegistryActor))
                    }
                  }
                )
              }
            )
          }
        )
      }
    }
}

object WesRunRoutes {
  def extractAuthorizationHeader: HttpHeader => Option[HttpHeader] = {
    case h: HttpHeader if h.name() == "Authorization" => Option(h)
    case _ => None
  }

  def completeCromwellResponse(future: => Future[WesResponse]): Route = {

    import WesResponseJsonSupport.WesResponseErrorFormat
    import cromwell.webservice.routes.wes.WesResponseJsonSupport.WesResponseFormat
    import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._

    onComplete(future) {
      case Success(response: WesResponse) => complete(response)
      case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
    }
  }

}

