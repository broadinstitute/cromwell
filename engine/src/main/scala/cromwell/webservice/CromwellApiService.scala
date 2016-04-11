package cromwell.webservice

import akka.actor._
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.webservice.WorkflowJsonSupport._
import lenthall.config.ScalaConfig._
import lenthall.spray.SwaggerUiResourceHttpService
import lenthall.spray.WrappedRoute._
import spray.http.MediaTypes._
import spray.http.StatusCodes
import spray.json._
import spray.routing.Directive.pimpApply
import spray.routing._

import scala.util.{Failure, Success, Try}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = "2.1.1"
}

object CromwellApiServiceActor {
  def props(workflowManagerActorRef: ActorRef, config: Config): Props = {
    Props(classOf[CromwellApiServiceActor], workflowManagerActorRef, config)
  }
}

class CromwellApiServiceActor(val workflowManager: ActorRef, config: Config)
  extends Actor with CromwellApiService with SwaggerService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  val possibleRoutes = workflowRoutes.wrapped("api", config.getBooleanOr("api.routeUnwrapped")) ~ swaggerUiResourceRoute

  def receive = runRoute(possibleRoutes)
}

trait CromwellApiService extends HttpService with PerRequestCreator {

  import CromwellApiServiceActor._

  val workflowManager: ActorRef

  private def invalidWorkflowId(id: String) = respondWithMediaType(`application/json`) {
    complete(StatusCodes.BadRequest, APIResponse.fail(new Throwable(s"Invalid workflow ID: '$id'.")).toJson.prettyPrint)
  }

  val workflowRoutes = queryRoute ~ workflowOutputsRoute ~ submitRoute ~ workflowStdoutStderrRoute ~ abortRoute ~
    callOutputsRoute ~ callStdoutStderrRoute ~ metadataRoute ~ timingRoute ~ callCachingRoute ~ statusRoute

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, workflowId) =>
      get {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowStatus(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def queryRoute =
    path("workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        get {
          requestContext =>
            perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowQuery(parameters))
        }
      }
    }

  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, workflowId) =>
      post {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowAbort(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def submitRoute =
    path("workflows" / Segment) { version =>
      post {
        formFields("wdlSource", "workflowInputs".?, "workflowOptions".?) { (wdlSource, workflowInputs, workflowOptions) =>
          requestContext =>
            val workflowSourceFiles = WorkflowSourceFiles(wdlSource, workflowInputs.getOrElse("{}"), workflowOptions.getOrElse("{}"))
            perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowSubmit(workflowSourceFiles))
        }
      }
    }

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, workflowId) =>
      get {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowOutputs(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def callOutputsRoute =
    path("workflows" / Segment / Segment / "outputs" / Segment) { (version, workflowId, callFqn) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerCallOutputs(w, callFqn))
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def callStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs" / Segment) { (version, workflowId, callFqn) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerCallStdoutStderr(w, callFqn))
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def workflowStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs") { (version, workflowId) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowStdoutStderr(w))
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def metadataRoute =
    path("workflows" / Segment / Segment / "metadata") { (version, workflowId) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowMetadata(w))
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def timingRoute =
    path("workflows" / Segment / Segment / "timing") { (version, workflowId) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(_) => getFromResource("workflowTimings/workflowTimings.html")
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def callCachingRoute =
    path("workflows" / Segment / Segment / "call-caching" ~ (Slash ~ Segment).?) { (version, workflowId, callFqn) =>
      parameterSeq { parameters =>
        val queryParameters = parameters map { case (k, v) => QueryParameter(k, v) }
        post {
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) =>
              requestContext =>
                perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerCallCaching(w, queryParameters, callFqn))
            case Failure(_) => invalidWorkflowId(workflowId)
          }
        }
      }
    }
}