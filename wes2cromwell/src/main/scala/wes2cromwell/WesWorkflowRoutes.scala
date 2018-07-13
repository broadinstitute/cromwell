package wes2cromwell

import akka.actor.{ActorRef, ActorSystem}
import akka.event.{Logging, LoggingAdapter}
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.{Directive1, Route}
import akka.http.scaladsl.server.directives.MethodDirectives.{delete, get, post}
import akka.http.scaladsl.server.directives.PathDirectives.path
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromiam.cromwell.CromwellClient
import wes2cromwell.WorkflowActor._
import net.ceedubs.ficus.Ficus._
import spray.json.JsObject
import cromiam.webservice.RequestSupport

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success}

// WorkflowRoutes implements the 'workflows' endpoint in WES
trait WesWorkflowRoutes extends JsonSupport with RequestSupport {
  // we leave these abstract, since they will be provided by the App
  implicit def system: ActorSystem

  val log: LoggingAdapter

  // other dependencies that Routes use
  def workflowActor: ActorRef

  def cromwellClient: CromwellClient

  // Make this configurable
  implicit val duration: FiniteDuration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit lazy val timeout: Timeout = duration

  lazy val workflowRoutes: Route =
    // TODO: factor the top of this into a path prefix in WesServer
    pathPrefix("ga4gh" / "wes" / "v1" / "workflows") {
      concat(
        pathEnd {
          concat(
            get {
              val futureWes: Future[Any] = workflowActor.ask(GetWorkflows)
              handleWesResponse(futureWes)
            },
            post {
              extractStrictRequest { request =>
                extractSubmission() { submission =>
                  complete { cromwellClient.forwardToCromwell(request.withEntity(submission.entity)) }
                }
              }
            }
          )
        },
        // workflows/{workflow_id}
        path(Segment) { workflowId =>
          concat(
            get {
              val maybeWorkflow: Future[Any] = workflowActor.ask(GetWorkflow(workflowId))
              handleWesResponse(maybeWorkflow)
            },
            delete {
              val futureWes: Future[Any] = workflowActor.ask(DeleteWorkflow(workflowId))
              handleWesResponse(futureWes)
            }
          )
        },
        // workflows/{workflow_id}/status
        path(Segment / "status") { workflowId =>
          get {
            val futureWes: Future[Any] = workflowActor.ask(GetWorkflowStatus(workflowId))
            handleWesResponse(futureWes)
          }
        }
      )
    }

  def extractSubmission(): Directive1[WesSubmission] = {
   formFields(
      "workflow_params",
      "workflow_type",
      "workflow_type_version",
      "tags".?,
      "workflow_engine_parameters".?,
      "workflow_url",
      "workflow_attachment".*
    ).as(WesSubmission)
  }

  // Common handler for some Wes Responses
  // TODO: understand if I can avoid re-constructing the responses
  def handleWesResponse(futureWes: Future[Any]): Route = {
    onComplete(futureWes.mapTo[WesResponse]) {
      case Success(wesResponse) => {
        wesResponse match {
          case WesResponseCreateWorkflowId(workflow_id) =>
            complete((StatusCodes.Created, WesResponseCreateWorkflowId(workflow_id)))
          case WesResponseDeleteWorkflowId(workflow_id) =>
            complete((StatusCodes.OK, WesResponseDeleteWorkflowId(workflow_id)))
          case WesResponseStatus(workflow_id, state) =>
            complete((StatusCodes.OK, WesResponseStatus(workflow_id, state)))
          case WesResponseWorkflowList(list) =>
            complete((StatusCodes.OK, WesResponseWorkflowList(list)))
          case WesResponseError(msg, status_code) =>
            complete((status_code, WesResponseError(msg, status_code)))
          case WesResponseWorkflowMetadata(workflowLog) =>
            complete((StatusCodes.OK, workflowLog))
        }
      }
      case Failure(ex) => {
        complete((
          StatusCodes.InternalServerError,
          WesResponseError(s"PostWorkflow exception: ${ex.getMessage}", StatusCodes.InternalServerError.intValue)
        ))
      }
    }
  }
}
