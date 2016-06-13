package cromwell.webservice

import akka.actor._
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryClient
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataBuilderActor
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
  def props(workflowManagerActorRef: ActorRef, config: Config): Props = {
    Props(new CromwellApiServiceActor(workflowManagerActorRef, config))
  }
}

class CromwellApiServiceActor(val workflowManager: ActorRef, config: Config)
  extends Actor with CromwellApiService with SwaggerService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  val possibleRoutes = workflowRoutes.wrapped("api", config.getBooleanOr("api.routeUnwrapped")) ~ swaggerUiResourceRoute

  def receive = runRoute(possibleRoutes)
}

trait CromwellApiService extends HttpService with PerRequestCreator with ServiceRegistryClient {
  val workflowManager: ActorRef

  private def invalidWorkflowId(id: String) = respondWithMediaType(`application/json`) {
    complete(StatusCodes.BadRequest, APIResponse.fail(new Throwable(s"Invalid workflow ID: '$id'.")).toJson.prettyPrint)
  }

  private def unsupportedEndpoint(endpoint: String) = respondWithMediaType(`application/json`) {
    complete(StatusCodes.BadRequest, APIResponse.fail(new UnsupportedOperationException(s"Unsupported endpoint: '$endpoint'. If this is unexpected, is your API version correct?")).toJson.prettyPrint)
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowStdoutStderrRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~
    callCachingRoute ~ statusRoute

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, workflowId) =>
      get {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, MetadataBuilderActor.props(serviceRegistryActor), GetStatus(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def queryRoute =
    path("workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        get {
          requestContext =>
            perRequest(requestContext, MetadataBuilderActor.props(serviceRegistryActor), WorkflowQuery(requestContext.request.uri, parameters))
        }
      }
    }

  def queryPostRoute =
    path("workflows" / Segment / "query") { version =>
      entity(as[Seq[Map[String, String]]]) { parameterMap =>
        post {
          requestContext =>
            perRequest(requestContext, MetadataBuilderActor.props(serviceRegistryActor), WorkflowQuery(requestContext.request.uri, parameterMap.flatMap(_.toSeq)))
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

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, workflowId) =>
      get {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            requestContext => perRequest(requestContext, MetadataBuilderActor.props(serviceRegistryActor), WorkflowOutputs(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
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
      parameters('outputs ? true, 'timings ? true).as(WorkflowMetadataQueryParameters) { parameters =>
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) =>
            version match {
              case _ => requestContext => perRequest(requestContext, MetadataBuilderActor.props(serviceRegistryActor), GetSingleWorkflowMetadataAction(w))
            }
          case Failure(_) => invalidWorkflowId(workflowId)
        }
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