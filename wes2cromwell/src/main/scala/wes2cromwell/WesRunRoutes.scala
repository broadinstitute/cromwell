package wes2cromwell

import java.net.URL

import akka.actor.{ActorRef, ActorSystem}
import akka.event.LoggingAdapter
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.{Directive1, Route}
import akka.http.scaladsl.server.directives.MethodDirectives.post
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import cromiam.webservice.RequestSupport
import wes2cromwell.WesResponseJsonSupport._
// import akka.pattern.ask

//import scala.concurrent.Future
import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global

import scala.util.{Failure, Success}

// WorkflowRoutes implements the 'runs' endpoint in WES
trait WesRunRoutes extends JsonSupport with RequestSupport {
  // we leave these abstract, since they will be provided by the App
  implicit def system: ActorSystem

  val log: LoggingAdapter

  // other dependencies that Routes use
  def workflowActor: ActorRef

  def cromwellUrl: URL
  def cromwellApiVersion: String
  def cromwellPath: URL = new URL(cromwellUrl.toString + s"/api/workflows/$cromwellApiVersion")

  // Make this configurable
  implicit val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

  lazy val wes2CromwellInterface = new Wes2CromwellInterface(cromwellPath)

  lazy val runRoutes: Route =
    // TODO: factor the top of this into a path prefix in WesServer
    pathPrefix("ga4gh" / "wes" / "v1" / "runs") {
      concat(
        pathEnd {
          concat(
            get {
              //val futureWes: Future[Any] = workflowActor.ask(GetWorkflows)
              //handleWesResponse(futureWes)
              complete("HI")
            },
            post {
              extractStrictRequest { request =>
                extractSubmission() { submission =>
                  onComplete(wes2CromwellInterface.submit(submission, request)) {
                    case Success(a) => complete(a)
                    case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
                  }
                }
              }
            }
          )
        },
        path(Segment) { workflowId =>
          concat(
            get {
              //val maybeWorkflow: Future[Any] = workflowActor.ask(GetWorkflow(workflowId))
              //handleWesResponse(maybeWorkflow)
              complete("HI")
            },
            delete {
              extractRequest { request =>
                onComplete(wes2CromwellInterface.abort(workflowId, request)) {
                  case Success(a) => complete(a)
                  case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
                }
              }
            }
          )
        },
        path(Segment / "status") { workflowId =>
          get {
            extractRequest { request =>
              onComplete(wes2CromwellInterface.status(workflowId, request)) {
                case Success(a) => complete(a)
                case Failure(e) => complete(WesErrorResponse(e.getMessage, StatusCodes.InternalServerError.intValue))
              }
            }
          }
        }
      )
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

  // Common handler for some Wes Responses
  // TODO: understand if I can avoid re-constructing the responses
//  def handleWesResponse(futureWes: Future[Any]): Route = {
//    onComplete(futureWes.mapTo[WesResponse]) {
//      case Success(wesResponse) => {
//        wesResponse match {
//          case WesResponseCreateWorkflowId(workflow_id) =>
//            complete((StatusCodes.Created, WesResponseCreateWorkflowId(workflow_id)))
//          case WesResponseDeleteWorkflowId(workflow_id) =>
//            complete((StatusCodes.OK, WesResponseDeleteWorkflowId(workflow_id)))
//          case WesResponseStatus(workflow_id, state) =>
//            complete((StatusCodes.OK, WesResponseStatus(workflow_id, state)))
//          case WesResponseWorkflowList(list) =>
//            complete((StatusCodes.OK, WesResponseWorkflowList(list)))
//          case WesResponseError(msg, status_code) =>
//            complete((status_code, WesResponseError(msg, status_code)))
//          case WesResponseWorkflowMetadata(workflowLog) =>
//            complete((StatusCodes.OK, workflowLog))
//        }
//      }
//      case Failure(ex) => {
//        complete((
//          StatusCodes.InternalServerError,
//          WesResponseError(s"PostWorkflow exception: ${ex.getMessage}", StatusCodes.InternalServerError.intValue)
//        ))
//      }
//    }
 // }
}
