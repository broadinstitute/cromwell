package cromwell.webservice

import akka.actor._
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.instrumentation.Instrumentation.Monitor
import cromwell.webservice.CromwellApiServiceActor.traceName
import cromwell.webservice.WorkflowJsonSupport._
import lenthall.config.ScalaConfig._
import lenthall.spray.SwaggerUiResourceHttpService
import lenthall.spray.WrappedRoute._
import spray.http.MediaTypes._
import spray.http.StatusCodes
import spray.json._
import spray.routing.Directive.pimpApply
import spray.routing._
import spray.httpx.SprayJsonSupport._

import scala.util.{Failure, Success, Try}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = "2.1.1"
}

object CromwellApiServiceActor {
  def props(workflowManagerActorRef: ActorRef, validateActorRef: ActorRef, config: Config): Props = {
    Props(classOf[CromwellApiServiceActor], workflowManagerActorRef, validateActorRef, config)
  }

  def traceName(name: String) = Monitor.traceName(name)
}

class CromwellApiServiceActor(val workflowManager: ActorRef, val workflowDescriptorMaterializer: ActorRef, config: Config)
  extends Actor with CromwellApiService with SwaggerService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  private val swaggerRoutes = traceName("swaggerResource") { swaggerResourceRoute } ~ traceName("swaggerUi") { swaggerUiRoute }
  val possibleRoutes = workflowRoutes.wrapped("api", config.getBooleanOr("api.routeUnwrapped")) ~ swaggerRoutes

  def receive = runRoute(possibleRoutes)
}

trait CromwellApiService extends HttpService with PerRequestCreator {
  val workflowManager: ActorRef
  val workflowDescriptorMaterializer: ActorRef

  private def invalidWorkflowId(id: String) = respondWithMediaType(`application/json`) {
    complete(StatusCodes.BadRequest, APIResponse.fail(new Throwable(s"Invalid workflow ID: '$id'.")).toJson.prettyPrint)
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowStdoutStderrRoute ~ abortRoute ~ callOutputsRoute ~ callStdoutStderrRoute ~ validateRoute ~ metadataRoute ~
    timingRoute ~ callCachingRoute ~ statusRoute

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, workflowId) =>
      traceName("status") {
        get {
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowStatus(w))
            case Failure(ex) => invalidWorkflowId(workflowId)
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
              perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowQuery(parameters))
          }
        }
      }
    }

  def queryPostRoute =
    path("workflows" / Segment / "query") { version =>
      entity(as[Seq[Map[String, String]]]) { parameterMap =>
        post {
          requestContext =>
            perRequest(requestContext, CromwellApiHandler.props(workflowManager),
              CromwellApiHandler.ApiHandlerWorkflowQuery(parameterMap.flatMap(_.toSeq)))
        }
      }
    }

  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, workflowId) =>
      traceName("abort") {
        post {
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowAbort(w))
            case Failure(ex) => invalidWorkflowId(workflowId)
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
              perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowSubmit(workflowSourceFiles))
          }
        }
      }
    }

  def submitBatchRoute =
    path("workflows" / Segment / "batch") { version =>
      post {
        formFields("wdlSource", "workflowInputs", "workflowOptions".?) {
          (wdlSource, workflowInputs, workflowOptions) =>
            requestContext =>
              import spray.json._
              workflowInputs.parseJson match {
                case JsArray(inputses) =>
                  val sources = inputses.map(inputs => WorkflowSourceFiles(wdlSource, inputs.compactPrint, workflowOptions.getOrElse("{}")))
                  perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowSubmitBatch(sources))
                case _ => reject
              }
        }
      }
    }

  def validateRoute =
    path("workflows" / Segment / "validate") { version =>
      traceName("validate") {
        post {
          formFields("wdlSource", "workflowInputs".?, "workflowOptions".?) { (wdlSource, workflowInputs, workflowOptions) =>
            requestContext =>
              perRequest(requestContext, CromwellApiHandler.props(workflowDescriptorMaterializer), CromwellApiHandler.ApiHandlerValidateWorkflow(WorkflowId.randomId(), wdlSource, workflowInputs, workflowOptions))
          }
        }
      }
    }

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, workflowId) =>
      traceName("workflowOutputs") {
        get {
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowOutputs(w))
            case Failure(ex) => invalidWorkflowId(workflowId)
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
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerCallOutputs(w, callFqn))
          case Failure(_) => invalidWorkflowId(workflowId)
        }
      }
    }

  def callStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs" / Segment) { (version, workflowId, callFqn) =>
      traceName("callLogs") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            // This currently does not attempt to parse the call name for conformation to any pattern.
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerCallStdoutStderr(w, callFqn))
          case Failure(_) => invalidWorkflowId(workflowId)
        }
      }
    }

  def workflowStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs") { (version, workflowId) =>
      traceName("workflowLogs") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.ApiHandlerWorkflowStdoutStderr(w))
          case Failure(_) => invalidWorkflowId(workflowId)
        }
      }
    }

  def metadataRoute =
    path("workflows" / Segment / Segment / "metadata") { (version, workflowId) =>
      traceName("workflowMetadata") {
        parameters('outputs ? true, 'timings ? true).as(WorkflowMetadataQueryParameters) { parameters =>
          Try(WorkflowId.fromString(workflowId)) match {
            case Success(w) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager),
                CromwellApiHandler.ApiHandlerWorkflowMetadata(w, parameters))
            case Failure(_) => invalidWorkflowId(workflowId)
          }
        }
      }
    }

  def timingRoute =
    path("workflows" / Segment / Segment / "timing") { (version, workflowId) =>
      traceName("timing") {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(_) => getFromResource("workflowTimings/workflowTimings.html")
          case Failure(_) => invalidWorkflowId(workflowId)
        }
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
