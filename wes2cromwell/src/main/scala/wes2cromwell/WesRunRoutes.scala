package wes2cromwell

import java.net.URL

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.model.{HttpHeader, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.{Directive1, Route}
import akka.http.scaladsl.server.directives.MethodDirectives.post
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import cromiam.webservice.RequestSupport
import wes2cromwell.WesResponseJsonSupport._
import WesRunRoutes._
import akka.stream.ActorMaterializer

import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Failure, Success}

trait WesRunRoutes extends RequestSupport {
  implicit def system: ActorSystem
  implicit def materializer: ActorMaterializer

  val log: LoggingAdapter

  def cromwellUrl: URL
  def cromwellApiVersion: String
  def cromwellPath: URL = new URL(cromwellUrl.toString + s"/api/workflows/$cromwellApiVersion")

  implicit lazy val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

  lazy val wes2CromwellInterface = new Wes2CromwellInterface(cromwellPath)

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
                      completeCromwellResponse(wes2CromwellInterface.listRuns(pageSize, pageToken, cromwellRequestHeaders))
                    }
                  },
                  post {
                    extractStrictRequest { request =>
                      extractSubmission() { submission =>
                        completeCromwellResponse(wes2CromwellInterface.runWorkflow(submission, cromwellRequestHeaders))
                      }
                    }
                  }
                )
              },
              path(Segment) { workflowId =>
                concat(
                  get {
                    completeCromwellResponse(wes2CromwellInterface.runLog(workflowId, cromwellRequestHeaders))
                  },
                  delete {
                    completeCromwellResponse(wes2CromwellInterface.cancelRun(workflowId, cromwellRequestHeaders))
                  }
                )
              },
              path(Segment / "status") { workflowId =>
                get {
                  completeCromwellResponse(wes2CromwellInterface.runStatus(workflowId, cromwellRequestHeaders))
                }
              }
            )
          }
        )
      }
    }

  def extractSubmission(): Directive1[WesSubmission] = {
    formFields((
      "workflow_params",
      "workflow_type",
      "workflow_type_version",
      "tags".?,
      "workflow_engine_parameters".?,
      "workflow_url",
      "workflow_attachment".as[String].*
    )).as(WesSubmission)
  }
}

object WesRunRoutes {
  def extractAuthorizationHeader: HttpHeader => Option[HttpHeader] = {
    case h: HttpHeader if h.name() == "Authorization" => Option(h)
    case _ => None
  }

  def completeCromwellResponse(future: ⇒ Future[WesResponse]): Route = {
    onComplete(future) {
      case Success(a) => complete(a)
      case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
    }
  }
}
