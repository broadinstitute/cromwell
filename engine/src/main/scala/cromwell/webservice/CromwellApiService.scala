package cromwell.webservice

import java.util.UUID

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.server.Directives._

import scala.concurrent.{ExecutionContext, Future}
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.Multipart.BodyPart
import akka.stream.ActorMaterializer
import cromwell.engine.backend.BackendConfiguration
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import cromwell.core.{WorkflowAborted, WorkflowId, WorkflowSubmitted}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreEngineActor, WorkflowStoreSubmitActor}
import akka.pattern.ask
import akka.util.{ByteString, Timeout}
import net.ceedubs.ficus.Ficus._
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.metadata.{MetadataBuilderActor, WorkflowQueryPagination}
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse, MetadataBuilderActorResponse}
import WorkflowJsonSupport._
import akka.http.scaladsl.coding.{Deflate, Gzip, NoCoding}
import akka.http.scaladsl.server.Route
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.labels.Labels
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{BuiltCallCacheDiffResponse, CachedCallNotFoundException, CallCacheDiffActorResponse, FailedCallCacheDiffResponse}
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheDiffActor, CallCacheDiffQueryParameter}
import cromwell.engine.workflow.workflowstore.WorkflowStoreEngineActor.WorkflowStoreEngineAbortResponse
import cromwell.webservice.LabelsManagerActor._
import lenthall.exception.AggregatedMessageException
import spray.json.JsObject

import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

trait CromwellApiService {
  import cromwell.webservice.CromwellApiService._

  implicit def actorRefFactory: ActorRefFactory
  implicit val materializer: ActorMaterializer
  implicit val ec: ExecutionContext

  val workflowStoreActor: ActorRef
  val workflowManagerActor: ActorRef
  val serviceRegistryActor: ActorRef

  // Derive timeouts (implicit and not) from akka http's request timeout since there's no point in being higher than that
  implicit val duration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit val timeout: Timeout = duration

  val routes =
    path("workflows" / Segment / "backends") { version =>
      get { complete(backendResponse) }
    } ~
    path("engine" / Segment / "stats") { version =>
      get {
        onComplete(workflowManagerActor.ask(WorkflowManagerActor.EngineStatsCommand).mapTo[EngineStatsActor.EngineStats]) {
          case Success(stats) => complete(stats)
          case Failure(_) => new RuntimeException("Unable to gather engine stats").failRequest(StatusCodes.InternalServerError)
        }
      }
    } ~
    path("engine" / Segment / "version") { version =>
      get { complete(versionResponse) }
    } ~
    path("workflows" / Segment / Segment / "status") { (version, possibleWorkflowId) =>
      get { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetStatus(w)) }
    } ~
    path("workflows" / Segment / Segment / "outputs") { (version, possibleWorkflowId) =>
      get { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => WorkflowOutputs(w)) }
    } ~
    path("workflows" / Segment / Segment / "logs") { (version, possibleWorkflowId) =>
      get { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetLogs(w)) }
    } ~
    path("workflows" / Segment / "query") { version =>
      (post | get) {
        parameterSeq { parameters =>
          extractUri { uri =>
            metadataQueryRequest(parameters, uri)
          }
        }
      }
    } ~
    encodeResponseWith(Gzip, Deflate, NoCoding) {
      path("workflows" / Segment / Segment / "metadata") { (version, possibleWorkflowId) =>
        parameters(('includeKey.*, 'excludeKey.*, 'expandSubWorkflows.as[Boolean].?)) { (includeKeys, excludeKeys, expandSubWorkflowsOption) =>
          val includeKeysOption = NonEmptyList.fromList(includeKeys.toList)
          val excludeKeysOption = NonEmptyList.fromList(excludeKeys.toList)
          val expandSubWorkflows = expandSubWorkflowsOption.getOrElse(false)

          (includeKeysOption, excludeKeysOption) match {
            case (Some(_), Some(_)) =>
              val e = new IllegalArgumentException("includeKey and excludeKey may not be specified together")
              e.failRequest(StatusCodes.BadRequest)
            case (_, _) => metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetSingleWorkflowMetadataAction(w, includeKeysOption, excludeKeysOption, expandSubWorkflows))
          }
        }
      }
    } ~
    path("workflows" / Segment / "callcaching" / "diff") { version =>
      parameterSeq { parameters =>
        get {
          CallCacheDiffQueryParameter.fromParameters(parameters) match {
            case Valid(queryParameter) =>
              val diffActor = actorRefFactory.actorOf(CallCacheDiffActor.props(serviceRegistryActor), "CallCacheDiffActor-" + UUID.randomUUID())
              onComplete(diffActor.ask(queryParameter).mapTo[CallCacheDiffActorResponse]) {
                case Success(r: BuiltCallCacheDiffResponse) => complete(r.response)
                case Success(r: FailedCallCacheDiffResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
                case Failure(e: CachedCallNotFoundException) => e.errorRequest(StatusCodes.NotFound)
                case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
              }
            case Invalid(errors) =>
              val e = AggregatedMessageException("Wrong parameters for call cache diff query", errors.toList)
              e.errorRequest(StatusCodes.BadRequest)
          }
        }
      }
    } ~
    path("workflows" / Segment / Segment / "timing") { (version, possibleWorkflowId) =>
      onComplete(validateWorkflowId(possibleWorkflowId)) {
        case Success(_) => getFromResource("workflowTimings/workflowTimings.html")
        case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
      }
    } ~
    path("workflows" / Segment / Segment / "abort") { (version, possibleWorkflowId) =>
      post {
        val response = validateWorkflowId(possibleWorkflowId) flatMap { w =>
          workflowStoreActor.ask(WorkflowStoreActor.AbortWorkflow(w, workflowManagerActor)).mapTo[WorkflowStoreEngineAbortResponse]
        }

        onComplete(response) {
          case Success(WorkflowStoreEngineActor.WorkflowAborted(id)) => complete(WorkflowAbortResponse(id.toString, WorkflowAborted.toString))
          case Success(WorkflowStoreEngineActor.WorkflowAbortFailed(_, e: IllegalStateException)) => e.errorRequest(StatusCodes.Forbidden)
          case Success(WorkflowStoreEngineActor.WorkflowAbortFailed(_, e: WorkflowNotFoundException)) => e.errorRequest(StatusCodes.NotFound)
          case Success(WorkflowStoreEngineActor.WorkflowAbortFailed(_, e)) => e.errorRequest(StatusCodes.InternalServerError)
          case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
          case Failure(e: InvalidWorkflowException) => e.failRequest(StatusCodes.BadRequest)
          case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
        }
      }
    } ~
    path("workflows" / Segment / Segment / "labels") { (version, possibleWorkflowId) =>
      entity(as[Map[String, String]]) { parameterMap =>
        patch {
          Labels.validateMapOfLabels(parameterMap) match {
            case Valid(labels) =>
              val response = validateWorkflowId(possibleWorkflowId) flatMap { id =>
                val lma = actorRefFactory.actorOf(LabelsManagerActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher))
                lma.ask(LabelsAddition(LabelsData(id, labels))).mapTo[LabelsManagerActorResponse]
              }
              onComplete(response) {
                case Success(r: BuiltLabelsManagerResponse) => complete(r.response)
                case Success(e: FailedLabelsManagerResponse) => e.reason.failRequest(StatusCodes.InternalServerError)
                case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)

              }
            case Invalid(e) =>
              val iae = new IllegalArgumentException(e.toList.mkString(","))
              iae.failRequest(StatusCodes.BadRequest)
          }
        }
      }
    } ~
  path("workflows" / Segment) { version =>
      post {
        entity(as[Multipart.FormData]) { formData =>
          submitRequest(formData, true)
        }
      }
    } ~
  path("workflows" / Segment / "batch") { version =>
    post {
      entity(as[Multipart.FormData]) { formData =>
        submitRequest(formData, false)
      }
    }
  }

  private def submitRequest(formData: Multipart.FormData, isSingleSubmission: Boolean): Route = {
    val allParts: Future[Map[String, ByteString]] = formData.parts.mapAsync[(String, ByteString)](1) {
      case b: BodyPart => b.toStrict(duration).map(strict => b.name -> strict.entity.data)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)

    onComplete(allParts) {
      case Success(data) =>
        PartialWorkflowSources.fromSubmitRoute(data, allowNoInputs = isSingleSubmission) match {
          case Success(workflowSourceFiles) if isSingleSubmission && workflowSourceFiles.size == 1 =>
            onComplete(workflowStoreActor.ask(WorkflowStoreActor.SubmitWorkflow(workflowSourceFiles.head)).mapTo[WorkflowStoreSubmitActor.WorkflowSubmittedToStore]) {
              case Success(w) => complete((StatusCodes.Created, WorkflowSubmitResponse(w.workflowId.toString, WorkflowSubmitted.toString)))
              case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
            }
          // Catches the case where someone has gone through the single submission endpoint w/ more than one workflow
          case Success(workflowSourceFiles) if isSingleSubmission =>
            val e = new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint.")
            e.failRequest(StatusCodes.BadRequest)
          case Success(workflowSourceFiles) =>
            onComplete(workflowStoreActor.ask(WorkflowStoreActor.BatchSubmitWorkflows(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList))).mapTo[WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore]) {
              case Success(w) =>
                val responses = w.workflowIds map { id => WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString) }
                complete((StatusCodes.Created, responses.toList))
              case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
            }
          case Failure(t) => t.failRequest(StatusCodes.BadRequest)
        }
      case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
    }
  }

  private def validateWorkflowId(possibleWorkflowId: String): Future[WorkflowId] = {
    Try(WorkflowId.fromString(possibleWorkflowId)) match {
      case Success(w) =>
        serviceRegistryActor.ask(ValidateWorkflowId(w)).mapTo[WorkflowValidationResponse] map {
          case RecognizedWorkflowId => w
          case UnrecognizedWorkflowId => throw UnrecognizedWorkflowException(s"Unrecognized workflow ID: $w")
          case FailedToCheckWorkflowId(t) => throw t
        }
      case Failure(t) => Future.failed(InvalidWorkflowException(s"Invalid workflow ID: '$possibleWorkflowId'."))
    }
  }

  private def metadataBuilderRequest(possibleWorkflowId: String, request: WorkflowId => ReadAction): Route = {
    val metadataBuilderActor = actorRefFactory.actorOf(MetadataBuilderActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName)
    val response = validateWorkflowId(possibleWorkflowId) flatMap { w => metadataBuilderActor.ask(request(w)).mapTo[MetadataBuilderActorResponse] }

    onComplete(response) {
      case Success(r: BuiltMetadataResponse) => complete(r.response)
      case Success(r: FailedMetadataResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
      case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
      case Failure(e: InvalidWorkflowException) => e.failRequest(StatusCodes.BadRequest)
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }

  protected[this] def metadataQueryRequest(parameters: Seq[(String, String)], uri: Uri): Route = {
    val response = serviceRegistryActor.ask(WorkflowQuery(parameters)).mapTo[MetadataQueryResponse]

    onComplete(response) {
      case Success(w: WorkflowQuerySuccess) =>
        val headers = WorkflowQueryPagination.generateLinkHeaders(uri, w.meta)
        respondWithHeaders(headers) {
          complete(w.response)
        }
      case Success(w: WorkflowQueryFailure) => w.reason.failRequest(StatusCodes.BadRequest)
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }
}

object CromwellApiService {
  import spray.json._

  implicit class EnhancedThrowable(val e: Throwable) extends AnyVal {
    def failRequest(statusCode: StatusCode): Route = complete((statusCode, APIResponse.fail(e).toJson.prettyPrint))
    def errorRequest(statusCode: StatusCode): Route = complete((statusCode, APIResponse.error(e).toJson.prettyPrint))
  }

  final case class BackendResponse(supportedBackends: List[String], defaultBackend: String)

  final case class UnrecognizedWorkflowException(message: String) extends Exception(message)
  final case class InvalidWorkflowException(message: String) extends Exception(message)

  val backendResponse = BackendResponse(BackendConfiguration.AllBackendEntries.map(_.name).sorted, BackendConfiguration.DefaultBackendEntry.name)
  val versionResponse = JsObject(Map("cromwell" -> ConfigFactory.load("cromwell-version.conf").getConfig("version").getString("cromwell").toJson))
}
