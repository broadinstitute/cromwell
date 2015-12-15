package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.engine.workflow.{ValidateActor, WorkflowOptions}
import cromwell.engine.{WorkflowId, WorkflowSourceFiles}
import cromwell.instrumentation.Instrumentation.Monitor
import cromwell.webservice.CromwellApiServiceActor.traceName
import lenthall.config.ScalaConfig._
import lenthall.spray.SwaggerUiResourceHttpService
import lenthall.spray.WrappedRoute._
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

  def traceName(name: String) = Monitor.traceName(name)
}

class CromwellApiServiceActor(val workflowManager: ActorRef, config: Config)
  extends Actor with CromwellApiService with SwaggerService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  private val swaggerRoutes = traceName("swaggerResource") { swaggerResourceRoute } ~ traceName("swaggerUi") { swaggerUiRoute }
  val possibleRoutes = workflowRoutes.wrapped("api", config.getBooleanOr("api.routeUnwrapped")) ~ swaggerRoutes

  def receive = runRoute(possibleRoutes)
}

trait CromwellApiService extends HttpService with PerRequestCreator {

  import CromwellApiServiceActor._

  val workflowManager: ActorRef

  val workflowRoutes = queryRoute ~ workflowOutputsRoute ~ submitRoute ~ workflowStdoutStderrRoute ~ abortRoute ~
    callOutputsRoute ~ callStdoutStderrRoute ~ validateRoute ~ metadataRoute ~ timingRoute ~ callCachingRoute

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, id) =>
      traceName("status") {
        get {
          Try(WorkflowId.fromString(id)) match {
            case Success(workflowId) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStatus(workflowId))
            case Failure(ex) => complete(StatusCodes.BadRequest)
          }
        }
      }
    }

  def queryRoute =
    path("workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        traceName("query") {
          get {
            requestContext =>
              perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowQuery(parameters))
          }
        }
      }
    }

  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, id) =>
      traceName("abort") {
        post {
          Try(WorkflowId.fromString(id)) match {
            case Success(workflowId) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowAbort(workflowId))
            case Failure(ex) => complete(StatusCodes.BadRequest)
          }
        }
      }
    }

  def submitRoute =
    path("workflows" / Segment) { version =>
      traceName("submit") {
        post {
          formFields("wdlSource", "workflowInputs".?, "workflowOptions".?) { (wdlSource, workflowInputs, workflowOptions) =>
            requestContext =>
              val workflowSourceFiles = WorkflowSourceFiles(wdlSource, workflowInputs.getOrElse("{}"), workflowOptions.getOrElse("{}"))
              perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowSubmit(workflowSourceFiles))
          }
        }
      }
    }

  def validateRoute =
    path("workflows" / Segment / "validate") { version =>
      traceName("validate") {
        post {
          formFields("wdlSource", "workflowInputs") { (wdlSource, workflowInputs) =>
            requestContext =>
              perRequest(
                requestContext,
                ValidateActor.props(wdlSource, workflowInputs),
                ValidateActor.ValidateWorkflow)
          }
        }
      }
    }

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, id) =>
      traceName("workflowOutputs") {
        get {
          Try(WorkflowId.fromString(id)) match {
            case Success(workflowId) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowOutputs(workflowId))
            case Failure(ex) => complete(StatusCodes.BadRequest)
          }
        }
      }
    }

  def callOutputsRoute =
    path("workflows" / Segment / Segment / "outputs" / Segment) { (version, workflowId, callFqn) =>
      traceName("callOutputs") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            // This currently does not attempt to parse the call name for conformation to any pattern.
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallOutputs(w, callFqn))
          case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
        }
      }
    }

  def callStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs" / Segment) { (version, workflowId, callFqn) =>
      traceName("callLogs") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            // This currently does not attempt to parse the call name for conformation to any pattern.
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallStdoutStderr(w, callFqn))
          case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
        }
      }
    }

  def workflowStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs") { (version, workflowId) =>
      traceName("workflowLogs") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStdoutStderr(w))
          case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
        }
      }
    }

  def metadataRoute =
    path("workflows" / Segment / Segment / "metadata") { (version, workflowId) =>
      traceName("workflowMetadata") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowMetadata(w))
          case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
        }
      }
    }

  def timingRoute =
    path("workflows" / Segment / Segment / "timing") { (version, workflowId) =>
      traceName("timing") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(_) => getFromResource("workflowTimings/workflowTimings.html")
          case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
        }
      }
    }

  def callCachingRoute =
    path("workflows" / Segment / Segment / "call-caching" ~ (Slash ~ Segment).?) { (version, workflowId, callFqn) =>
      parameterSeq { parameters =>
        val queryParameters = parameters map { case (k, v) => QueryParameter(k, v) }
        post { requestContext =>
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) => perRequest(requestContext, CromwellApiHandler.props(workflowManager),
              CromwellApiHandler.CallCaching(w, queryParameters, callFqn))
            case Failure(_) => complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
          }
        }
      }
    }
}
