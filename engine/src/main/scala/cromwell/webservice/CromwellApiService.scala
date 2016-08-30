package cromwell.webservice

import akka.actor._
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.engine.backend.BackendConfiguration
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataBuilderActor
import lenthall.spray.SwaggerUiResourceHttpService
import spray.http.MediaTypes._
import spray.http._
import spray.httpx.SprayJsonSupport._
import spray.json._
import spray.routing._

import scalaz.NonEmptyList

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = "2.1.1"
}

trait CromwellApiService extends HttpService with PerRequestCreator {
  val workflowManagerActor: ActorRef
  val workflowStoreActor: ActorRef
  val serviceRegistryActor: ActorRef

  def metadataBuilderProps: Props = MetadataBuilderActor.props(serviceRegistryActor)

  def handleMetadataRequest(message: AnyRef): Route = {
    requestContext =>
      perRequest(requestContext, metadataBuilderProps, message)
  }

  private def failBadRequest(exception: Exception, statusCode: StatusCode = StatusCodes.BadRequest) = respondWithMediaType(`application/json`) {
    complete(statusCode, APIResponse.fail(exception).toJson.prettyPrint)
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ callCachingRoute ~ statusRoute ~ backendRoute

  private def withRecognizedWorkflowId(possibleWorkflowId: String)(recognizedWorkflowId: WorkflowId => Route): Route = {
    def callback(requestContext: RequestContext) = new ValidationCallback {
      // The submitted value is malformed as a UUID and therefore not possibly recognized.
      override def onMalformed(possibleWorkflowId: String): Unit = {
        val exception = new RuntimeException(s"Invalid workflow ID: '$possibleWorkflowId'.")
        failBadRequest(exception)(requestContext)
      }

      override def onUnrecognized(possibleWorkflowId: String): Unit = {
        val exception = new RuntimeException(s"Unrecognized workflow ID: $possibleWorkflowId")
        failBadRequest(exception, StatusCodes.NotFound)(requestContext)
      }

      override def onFailure(possibleWorkflowId: String, throwable: Throwable): Unit = {
        val exception = new RuntimeException(s"Failed lookup attempt for workflow ID $possibleWorkflowId", throwable)
        failBadRequest(exception)(requestContext)
      }

      override def onRecognized(workflowId: WorkflowId): Unit = {
        recognizedWorkflowId(workflowId)(requestContext)
      }
    }

    requestContext => {
      val message = ValidateWorkflowIdAndExecute(possibleWorkflowId, callback(requestContext))
      serviceRegistryActor ! message
    }
  }

  def statusRoute =
    path("workflows" / Segment / Segment / "status") { (version, possibleWorkflowId) =>
      get {
        withRecognizedWorkflowId(possibleWorkflowId) { id =>
          handleMetadataRequest(GetStatus(id))
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
    path("workflows" / Segment / Segment / "abort") { (version, possibleWorkflowId) =>
      post {
        withRecognizedWorkflowId(possibleWorkflowId) { id =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManagerActor), CromwellApiHandler.ApiHandlerWorkflowAbort(id))
        }
      }
    }

  def submitRoute =
    path("workflows" / Segment) { version =>
      post {
        formFields("wdlSource", "workflowInputs".?, "workflowOptions".?) { (wdlSource, workflowInputs, workflowOptions) =>
          requestContext =>
            val workflowSourceFiles = WorkflowSourceFiles(wdlSource, workflowInputs.getOrElse("{}"), workflowOptions.getOrElse("{}"))
            perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowSubmit(workflowSourceFiles))
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
                case JsArray(Seq(x, xs@_*)) =>
                  val nelInputses = NonEmptyList.nel(x, xs.toList)
                  val sources = nelInputses.map(inputs => WorkflowSourceFiles(wdlSource, inputs.compactPrint, workflowOptions.getOrElse("{}")))
                  perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowSubmitBatch(sources))
                case JsArray(_) => failBadRequest(new RuntimeException("Nothing was submitted"))
                case _ => reject
              }
        }
      }
    }

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, possibleWorkflowId) =>
      get {
        withRecognizedWorkflowId(possibleWorkflowId) { id =>
          handleMetadataRequest(WorkflowOutputs(id))
        }
      }
    }

  def workflowLogsRoute =
    path("workflows" / Segment / Segment / "logs") { (version, possibleWorkflowId) =>
      get {
        withRecognizedWorkflowId(possibleWorkflowId) { id =>
          handleMetadataRequest(GetLogs(id))
        }
      }
    }

  def metadataRoute =
    path("workflows" / Segment / Segment / "metadata") { (version, possibleWorkflowId) =>
      parameterMultiMap { parameters =>
        // import scalaz_ & Scalaz._ add too many slow implicits, on top of the spray and json implicits
        import scalaz.syntax.std.list._
        val includeKeysOption = parameters.getOrElse("includeKey", List.empty).toNel
        val excludeKeysOption = parameters.getOrElse("excludeKey", List.empty).toNel
        (includeKeysOption, excludeKeysOption) match {
          case (Some(_), Some(_)) =>
            failBadRequest(new IllegalArgumentException("includeKey and excludeKey may not be specified together"))
          case _ =>
            withRecognizedWorkflowId(possibleWorkflowId) { id =>
              handleMetadataRequest(GetSingleWorkflowMetadataAction(id, includeKeysOption, excludeKeysOption))
            }
        }
      }
    }

  def timingRoute =
    path("workflows" / Segment / Segment / "timing") { (version, possibleWorkflowId) =>
      withRecognizedWorkflowId(possibleWorkflowId) { id =>
        getFromResource("workflowTimings/workflowTimings.html")
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

  def backendRoute =
    path("workflows" / Segment / "backends") { version =>
      get {
        complete {
          // Note that this is not using our standard per-request scheme, since the result is pre-calculated already
          backendResponse
        }
      }
    }

  val backendResponse = JsObject(Map(
    "supportedBackends" -> BackendConfiguration.AllBackendEntries.map(_.name).sorted.toJson,
    "defaultBackend" -> BackendConfiguration.DefaultBackendEntry.name.toJson
  ))

}

