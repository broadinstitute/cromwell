package cromwell.webservice

import akka.actor._
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.{WorkflowId, WorkflowOptionsJson, WorkflowSourceFilesCollection}
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
import wdl4s.{WdlJson, WdlSource}

import scala.util.{Failure, Success, Try}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = "2.1.1"
}

trait CromwellApiService extends HttpService with PerRequestCreator {
  val workflowManagerActor: ActorRef
  val workflowStoreActor: ActorRef
  val serviceRegistryActor: ActorRef

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

  private def failBadRequest(t: Throwable, statusCode: StatusCode = StatusCodes.BadRequest) = respondWithMediaType(`application/json`) {
    complete((statusCode, APIResponse.fail(t).toJson.prettyPrint))
  }

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ statsRoute ~ versionRoute

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
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowStoreActor), CromwellApiHandler.ApiHandlerWorkflowAbort(id, workflowManagerActor))
        }
      }
    }

  case class PartialWorkflowSources(wdlSource: Option[WdlSource], workflowInputs: Seq[WdlJson], workflowInputsAux: Map[Int, WdlJson], workflowOptions: Option[WorkflowOptionsJson], zippedImports: Option[Array[Byte]])
  object PartialWorkflowSources {
    private def workflowInputs(bodyPart: BodyPart): Seq[WdlJson] = {
      import spray.json._
      bodyPart.entity.data.asString.parseJson match {
        case JsArray(Seq(x, xs@_*)) => (List(x) ++ xs).map(_.compactPrint)
        case JsArray(_) => Seq.empty
        case v: JsValue => Seq(v.compactPrint)
      }
    }

    def partialSourcesToSourceCollections(partialSources: Try[PartialWorkflowSources], allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
      partialSources flatMap {
        case PartialWorkflowSources(Some(wdlSource), workflowInputs, workflowInputsAux, workflowOptions, wdlDependencies) =>
          //The order of addition allows for the expected override of colliding keys.
          val sortedInputAuxes = workflowInputsAux.toSeq.sortBy(_._1).map(x => Option(x._2))
          val wfInputs: Try[Seq[WdlJson]] = if (workflowInputs.isEmpty) {
            if (allowNoInputs) Success(Seq("{}")) else Failure(new IllegalArgumentException("No inputs were provided"))
          } else Success(workflowInputs map { workflowInputSet =>
            mergeMaps(Seq(Option(workflowInputSet)) ++ sortedInputAuxes).toString
          })
          wfInputs.map(_.map(x => WorkflowSourceFilesCollection(wdlSource, x, workflowOptions.getOrElse("{}"), wdlDependencies)))
        case other => Failure(new IllegalArgumentException(s"Incomplete workflow submission: $other"))
      }
    }

    def fromSubmitRoute(formData: MultipartFormData, allowNoInputs: Boolean): Try[Seq[WorkflowSourceFilesCollection]] = {
      val partialSources = Try(formData.fields.foldLeft(PartialWorkflowSources(None, Seq.empty, Map.empty, None, None)) { (partialSources: PartialWorkflowSources, bodyPart: BodyPart) =>
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
        } else {
          throw new IllegalArgumentException(s"Unexpected body part name: ${bodyPart.name.getOrElse("None")}")
        }
      })
      partialSourcesToSourceCollections(partialSources, allowNoInputs)
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

  def metadataRoute =
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

