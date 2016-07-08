package cromwell.webservice

import java.util.concurrent.TimeoutException

import akka.actor._
import com.typesafe.config.Config
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.services.MetadataServiceActor._
import cromwell.services.ServiceRegistryClient
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataBuilderActor
import lenthall.config.ScalaConfig._
import lenthall.spray.SwaggerUiResourceHttpService
import lenthall.spray.WrappedRoute._
import spray.http.HttpHeaders.`Content-Type`
import spray.http.MediaTypes._
import spray.http._
import spray.json._
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
  val timeoutError = APIResponse.error(new TimeoutException("The server was not able to produce a timely response to your request.")).toJson.prettyPrint

  def receive = handleTimeouts orElse runRoute(possibleRoutes)

  def handleTimeouts: Receive = {
    case Timedout(_: HttpRequest) =>
      sender() ! HttpResponse(StatusCodes.InternalServerError, HttpEntity(ContentType(MediaTypes.`application/json`), timeoutError)).withHeaders(`Content-Type`(`application/json`))
  }
}

trait CromwellApiService extends HttpService with PerRequestCreator with ServiceRegistryClient {
  val workflowManager: ActorRef
  def metadataBuilderProps: Props = MetadataBuilderActor.props(serviceRegistryActor)
  def handleMetadataRequest(message: AnyRef)(requestContext: RequestContext): Unit = {
    perRequest(requestContext, metadataBuilderProps, message)
  }

  private def invalidWorkflowId(id: String) = failBadRequest(new RuntimeException(s"Invalid workflow ID: '$id'."))

  private def failBadRequest(exception: Exception) = respondWithMediaType(`application/json`) {
    complete(StatusCodes.BadRequest, APIResponse.fail(exception).toJson.prettyPrint)
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowStdoutStderrRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~
    callCachingRoute ~ statusRoute

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, workflowId) =>
      get {
        Try(WorkflowId.fromString(workflowId)) match {
          case Success(w) => handleMetadataRequest(GetStatus(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def queryRoute =
    path("workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        get {
          requestContext =>
            perRequest(requestContext, metadataBuilderProps, WorkflowQuery(requestContext.request.uri, parameters))
        }
      }
    }

  def queryPostRoute =
    path("workflows" / Segment / "query") { version =>
      entity(as[Seq[Map[String, String]]]) { parameterMap =>
        post {
          requestContext =>
            perRequest(requestContext, metadataBuilderProps, WorkflowQuery(requestContext.request.uri, parameterMap.flatMap(_.toSeq)))
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
          case Success(w) => handleMetadataRequest(WorkflowOutputs(w))
          case Failure(ex) => invalidWorkflowId(workflowId)
        }
      }
    }

  def workflowStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs") { (version, workflowId) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) => handleMetadataRequest(GetLogs(w))
        case Failure(_) => invalidWorkflowId(workflowId)
      }
    }

  def metadataRoute =
    path("workflows" / Segment / Segment / "metadata") { (version, workflowId) =>
      parameterMultiMap { parameters =>
        // import scalaz_ & Scalaz._ add too many slow implicits, on top of the spray and json implicits
        import scalaz.syntax.std.list._
        val includeKeysOption = parameters.getOrElse("includeKey", List.empty).toNel
        val excludeKeysOption = parameters.getOrElse("excludeKey", List.empty).toNel
        (includeKeysOption, excludeKeysOption) match {
          case (Some(_), Some(_)) =>
            failBadRequest(new IllegalArgumentException("includeKey and excludeKey may not be specified together"))
          case _ =>
            Try(WorkflowId.fromString(workflowId)) match {
              case Success(workflowIdFromString) =>
                version match {
                  case _ => handleMetadataRequest(
                        GetSingleWorkflowMetadataAction(workflowIdFromString, includeKeysOption, excludeKeysOption))
                }
              case Failure(_) => invalidWorkflowId(workflowId)
            }
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
          // TODO: PBE: Certainly want to do something for this! But probably not to the WMA
          failBadRequest(new UnsupportedOperationException(s"Call caching is currently unsupported."))
        }
      }
    }
}