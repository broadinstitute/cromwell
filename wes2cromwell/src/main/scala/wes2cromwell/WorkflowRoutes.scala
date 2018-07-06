package wes2cromwell

import akka.actor.{ActorRef, ActorSystem}
import akka.event.Logging
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.server.directives.MethodDirectives.{delete, get, post}
import akka.http.scaladsl.server.directives.PathDirectives.path
import akka.http.scaladsl.server.directives.RouteDirectives.complete
import akka.pattern.ask
import akka.util.Timeout
import wes2cromwell.WorkflowActor._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success}

// WorkflowRoutes implements the 'workflows' endpoint in WES
trait WorkflowRoutes extends JsonSupport {
  // we leave these abstract, since they will be provided by the App
  implicit def system: ActorSystem

  lazy val log = Logging(system, classOf[WorkflowRoutes])

  // other dependencies that Routes use
  def workflowActor: ActorRef

  // Required by the `ask` (?) method below
  implicit lazy val timeout = Timeout(30.seconds) // usually we'd obtain the timeout from the system's configuration

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
              entity(as[WorkflowRequest]) { workflowRequest =>
                val futureWes: Future[Any] = workflowActor.ask(PostWorkflow(workflowRequest))
                handleWesResponse(futureWes)
              }
            }
          )
        },
        // workflows/{workflow_id}
        path(Segment) { workflowId =>
          concat(
            get {
              val maybeWorkflow: Future[Option[WorkflowLog]] =
                (workflowActor ? GetWorkflow(workflowId)).mapTo[Option[WorkflowLog]]
              rejectEmptyResponse {
                complete(maybeWorkflow)
              }
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

  // Common handler for some Wes Responses
  // TODO: understand if I can avoid re-constructing the responses
  def handleWesResponse(futureWes: Future[Any]) = {
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
            complete((status_code, WesResponseError(msg, status_code)) )
        }
      }
      case Failure(ex) => {
        complete((StatusCodes.InternalServerError, WesResponseError(s"PostWorkflow exception: ${ex.getMessage}", StatusCodes.InternalServerError.intValue)))
      }
    }
  }

}
