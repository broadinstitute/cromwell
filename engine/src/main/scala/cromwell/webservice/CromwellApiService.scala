package cromwell.webservice

import akka.actor._
import cats.data.Validated.{Invalid, Valid}
import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.cartesian._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowOptionsJson, WorkflowSourceFilesCollection}
import cromwell.engine.backend.BackendConfiguration
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataBuilderActor
import cromwell.core._
import cromwell.webservice.CromwellApiHandler.ApiHandlerCallCachingDiff
import lenthall.validation.ErrorOr.ErrorOr
import spray.http.MediaTypes._
import spray.http._
import spray.httpx.SprayJsonSupport._
import spray.json._
import spray.routing._
import wdl4s.{WdlJson, WdlSource}

import scala.util.{Failure, Success, Try}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = "2.1.1"
}

trait CromwellApiService extends HttpService with PerRequestCreator {
  def workflowManagerActor: ActorRef
  def workflowStoreActor: ActorRef
  def serviceRegistryActor: ActorRef
  def callCacheReadActor: ActorRef

  def toMap(someInput: Option[String]): Map[String, JsValue] = {
    import spray.json._
    someInput match {
      case Some(inputs: String) => inputs.parseJson match {
        case JsObject(inputMap) => inputMap
        case _ =>
          throw new RuntimeException(s"Submitted inputs couldn't be processed, please check for syntactical errors")
      }
      case None => Map.empty
    }
  }

  def mergeMaps(allInputs: Seq[Option[String]]): JsObject = {
    val convertToMap = allInputs.map(x => toMap(x))
    JsObject(convertToMap reduce (_ ++ _))
  }

  def metadataBuilderProps: Props = MetadataBuilderActor.props(serviceRegistryActor)

  def handleMetadataRequest(message: AnyRef): Route = {
    requestContext =>
      perRequest(requestContext, metadataBuilderProps, message)
  }

  def handleQueryMetadataRequest(parameters: Seq[(String, String)]): Route = {
    requestContext =>
      perRequest(requestContext, metadataBuilderProps, WorkflowQuery(requestContext.request.uri, parameters))
  }

  def handleCallCachingDiffRequest(parameters: Seq[(String, String)]): Route = {
    requestContext =>
      perRequest(requestContext, CromwellApiHandler.props(callCacheReadActor), ApiHandlerCallCachingDiff(parameters))
  }

  protected def failBadRequest(t: Throwable, statusCode: StatusCode = StatusCodes.BadRequest) = respondWithMediaType(`application/json`) {
    complete((statusCode, APIResponse.fail(t).toJson.prettyPrint))
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~ callCachingDiffRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ statsRoute ~ versionRoute

  protected def withRecognizedWorkflowId(possibleWorkflowId: String)(recognizedWorkflowId: WorkflowId => Route): Route = {
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
          handleQueryMetadataRequest(parameters)
        }
      }
    }

  def queryPostRoute =
    path("workflows" / Segment / "query") { version =>
      entity(as[Seq[Map[String, String]]]) { parameterMap =>
        post {
          handleQueryMetadataRequest(parameterMap.flatMap(_.toSeq))
        }
      }
    }

  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, possibleWorkflowId) =>
      post {
        withRecognizedWorkflowId(possibleWorkflowId) { id =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowAbort(id, workflowManagerActor))
        }
      }
    }

  def callCachingDiffRoute =
    path("workflows" / Segment / "callcaching") { version =>
      parameterSeq { parameters =>
        get {
          handleCallCachingDiffRequest(parameters)
        }
      }
    }

  case class PartialWorkflowSources
  (
    wdlSource: Option[WdlSource],
    workflowInputs: Vector[WdlJson],
    workflowInputsAux: Map[Int, WdlJson],
    workflowOptions: Option[WorkflowOptionsJson],
    customLabels: Option[WdlJson],
    zippedImports: Option[Array[Byte]])
  object PartialWorkflowSources {
    private def workflowInputs(bodyPart: BodyPart): Vector[WdlJson] = {
      import spray.json._
      bodyPart.entity.data.asString.parseJson match {
        case JsArray(Seq(x, xs@_*)) => (Vector(x) ++ xs).map(_.compactPrint)
        case JsArray(_) => Vector.empty
        case v: JsValue => Vector(v.compactPrint)
      }
    }

    def partialSourcesToSourceCollections(partialSources: ErrorOr[PartialWorkflowSources], allowNoInputs: Boolean): ErrorOr[Seq[WorkflowSourceFilesCollection]] = {

      def validateInputs(pws: PartialWorkflowSources): ErrorOr[Seq[WdlJson]] =
        (pws.workflowInputs.isEmpty, allowNoInputs) match {
          case (true, true) => Vector("{}").validNel
          case (true, false) => "No inputs were provided".invalidNel
          case _ =>
            val sortedInputAuxes = pws.workflowInputsAux.toSeq.sortBy { case (index, _) => index } map { case(_, inputJson) => Option(inputJson) }
            (pws.workflowInputs map { workflowInputSet: WdlJson => mergeMaps(Seq(Option(workflowInputSet)) ++ sortedInputAuxes).toString }).validNel
      }

      def validateOptions(options: Option[WorkflowOptionsJson]): ErrorOr[WorkflowOptions] =
        WorkflowOptions.fromJsonString(options.getOrElse("{}")).tryToErrorOr leftMap { _ map { i => s"Invalid workflow options provided: $i" } }

      def validateWdlSource(partialSource: PartialWorkflowSources): ErrorOr[WdlJson] = partialSource.wdlSource match {
        case Some(src) => src.validNel
        case _ => s"Incomplete workflow submission: $partialSource".invalidNel
      }

      partialSources match {
        case Valid(partialSource) =>
          (validateWdlSource(partialSource) |@| validateInputs(partialSource) |@| validateOptions(partialSource.workflowOptions)) map {
            case (wdlSource, wfInputs, wfOptions) =>
              wfInputs.map(x => WorkflowSourceFilesCollection(wdlSource, x, wfOptions.asPrettyJson, partialSource.customLabels.getOrElse("{}"), partialSource.zippedImports))
          }
        case Invalid(err) => err.invalid
      }
    }

    def fromSubmitRoute(formData: MultipartFormData, allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
      val partialSources = Try(formData.fields.foldLeft(PartialWorkflowSources(None, Vector.empty, Map.empty, None, None, None)) { (partialSources: PartialWorkflowSources, bodyPart: BodyPart) =>
        if (bodyPart.name.contains("wdlSource")) {
          partialSources.copy(wdlSource = Some(bodyPart.entity.data.asString))
        } else if (bodyPart.name.contains("workflowInputs")) {
          partialSources.copy(workflowInputs = workflowInputs(bodyPart))
        } else if (bodyPart.name.forall(_.startsWith("workflowInputs_"))) {
          val index = bodyPart.name.get.stripPrefix("workflowInputs_").toInt
          partialSources.copy(workflowInputsAux = partialSources.workflowInputsAux + (index -> bodyPart.entity.data.asString))
        } else if (bodyPart.name.contains("workflowOptions")) {
          partialSources.copy(workflowOptions = Some(bodyPart.entity.data.asString))
        } else if (bodyPart.name.contains("wdlDependencies")) {
          partialSources.copy(zippedImports = Some(bodyPart.entity.data.toByteArray))
        } else if (bodyPart.name.contains("customLabels")) {
          partialSources.copy(customLabels = Some(bodyPart.entity.data.asString))
        } else {
          throw new IllegalArgumentException(s"Unexpected body part name: ${bodyPart.name.getOrElse("None")}")
        }
      })
      partialSourcesToSourceCollections(partialSources.tryToErrorOr, allowNoInputs).errorOrToTry
    }
  }

  def submitRoute =
    path("workflows" / Segment) { version =>
      post {
        entity(as[MultipartFormData]) { formData =>
          PartialWorkflowSources.fromSubmitRoute(formData, allowNoInputs = true) match {
            case Success(workflowSourceFiles) if workflowSourceFiles.size == 1 =>
              requestContext => {
                perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowSubmit(workflowSourceFiles.head))
              }
            case Success(workflowSourceFiles) =>
              failBadRequest(new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint."))
            case Failure(t) =>
              failBadRequest(t)
          }
        }
      }
    }

  def submitBatchRoute =
    path("workflows" / Segment / "batch") { version =>
      post {
        entity(as[MultipartFormData]) { formData =>
          PartialWorkflowSources.fromSubmitRoute(formData, allowNoInputs = false) match {
            case Success(workflowSourceFiles) =>
              requestContext => {
                perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowSubmitBatch(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList)))
              }
            case Failure(t) =>
              failBadRequest(t)
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

  def metadataRoute = compressResponse() {
    path("workflows" / Segment / Segment / "metadata") { (version, possibleWorkflowId) =>
      parameterMultiMap { parameters =>
        val includeKeysOption = NonEmptyList.fromList(parameters.getOrElse("includeKey", List.empty))
        val excludeKeysOption = NonEmptyList.fromList(parameters.getOrElse("excludeKey", List.empty))
        val expandSubWorkflowsOption = {
          parameters.get("expandSubWorkflows") match {
            case Some(v :: Nil) => Try(v.toBoolean)
            case _ => Success(false)
          }
        }

        (includeKeysOption, excludeKeysOption, expandSubWorkflowsOption) match {
          case (Some(_), Some(_), _) =>
            failBadRequest(new IllegalArgumentException("includeKey and excludeKey may not be specified together"))
          case (_, _, Success(expandSubWorkflows)) =>
            withRecognizedWorkflowId(possibleWorkflowId) { id =>
              handleMetadataRequest(GetSingleWorkflowMetadataAction(id, includeKeysOption, excludeKeysOption, expandSubWorkflows))
            }
          case (_, _, Failure(ex)) => failBadRequest(new IllegalArgumentException(ex))
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

  def statsRoute =
    path("engine" / Segment / "stats") { version =>
      get {
        requestContext =>
          perRequest(requestContext, CromwellApiHandler.props(workflowManagerActor), CromwellApiHandler.ApiHandlerEngineStats)
      }
    }

  def versionRoute =
    path("engine" / Segment / "version") { version =>
      get {
        complete {
          lazy val versionConf = ConfigFactory.load("cromwell-version.conf").getConfig("version")
          versionResponse(versionConf)
        }
      }
    }

  def versionResponse(versionConf: Config) = JsObject(Map(
    "cromwell" -> versionConf.getString("cromwell").toJson
  ))

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

