package cromwell.webservice

import java.util.UUID

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.{ToEntityMarshaller, ToResponseMarshallable}
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.RawHeader
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.{AskTimeoutException, ask}
import akka.stream.ActorMaterializer
import akka.util.{ByteString, Timeout}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import common.exception.AggregatedMessageException
import common.util.VersionUtil
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.abort.{AbortResponse, WorkflowAbortFailureResponse, WorkflowAbortingResponse}
import cromwell.core.labels.Labels
import cromwell.core.{WorkflowAborting, WorkflowId, WorkflowSubmitted}
import cromwell.engine.backend.BackendConfiguration
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{AbortWorkflowCommand, WorkflowNotFoundException}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{BuiltCallCacheDiffResponse, CachedCallNotFoundException, CallCacheDiffActorResponse, FailedCallCacheDiffResponse}
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheDiffActor, CallCacheDiffQueryParameter}
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreSubmitActor}
import cromwell.server.CromwellShutdown
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{GetCurrentStatus, StatusCheckResponse}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.LabelsManagerActor._
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse, MetadataBuilderActorResponse}
import cromwell.webservice.metadata.{MetadataBuilderActor, WorkflowQueryPagination}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

trait CromwellApiService extends HttpInstrumentation {
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

  val engineRoutes = concat(
    path("engine" / Segment / "stats") { _ =>
      get {
        onComplete(workflowManagerActor.ask(WorkflowManagerActor.EngineStatsCommand).mapTo[EngineStatsActor.EngineStats]) {
          case Success(stats) => complete(ToResponseMarshallable(stats))
          case Failure(_) => new RuntimeException("Unable to gather engine stats").failRequest(StatusCodes.InternalServerError)
        }
      }
    },
    path("engine" / Segment / "version") { _ =>
      get { complete(versionResponse) }
    },
    path("engine" / Segment / "status") { _ =>
      onComplete(serviceRegistryActor.ask(GetCurrentStatus).mapTo[StatusCheckResponse]) {
        case Success(status) =>
          val httpCode = if (status.ok) StatusCodes.OK else StatusCodes.InternalServerError
          complete(ToResponseMarshallable((httpCode, status.systems)))
        case Failure(_) => new RuntimeException("Unable to gather engine status").failRequest(StatusCodes.InternalServerError)
      }
    }
  )

  val workflowRoutes =
    path("workflows" / Segment / "backends") { _ =>
      get { instrumentRequest { complete(ToResponseMarshallable(backendResponse)) } }
    } ~
    path("workflows" / Segment / Segment / "status") { (_, possibleWorkflowId) =>
      get {  instrumentRequest { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetStatus(w)) } }
    } ~
    path("workflows" / Segment / Segment / "outputs") { (_, possibleWorkflowId) =>
      get {  instrumentRequest { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => WorkflowOutputs(w)) } }
    } ~
    path("workflows" / Segment / Segment / "logs") { (_, possibleWorkflowId) =>
      get {  instrumentRequest { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetLogs(w)) } }
    } ~
    path("workflows" / Segment / "query") { _ =>
      get {
        instrumentRequest {
          parameterSeq { parameters =>
            extractUri { uri =>
              metadataQueryRequest(parameters, uri)
            }
          }
        }
      } ~
      post {
        instrumentRequest {
          entity(as[Seq[Map[String, String]]]) { parameterMap =>
            extractUri { uri =>
              metadataQueryRequest(parameterMap.flatMap(_.toSeq), uri)
            }
          }
        }
      }
    } ~
    encodeResponse {
      path("workflows" / Segment / Segment / "metadata") { (_, possibleWorkflowId) =>
        instrumentRequest {
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
      }
    } ~
    path("workflows" / Segment / "callcaching" / "diff") { _ =>
      parameterSeq { parameters =>
        get {
          instrumentRequest {
            CallCacheDiffQueryParameter.fromParameters(parameters) match {
              case Valid(queryParameter) =>
                val diffActor = actorRefFactory.actorOf(CallCacheDiffActor.props(serviceRegistryActor), "CallCacheDiffActor-" + UUID.randomUUID())
                onComplete(diffActor.ask(queryParameter).mapTo[CallCacheDiffActorResponse]) {
                  case Success(r: BuiltCallCacheDiffResponse) => complete(r.response)
                  case Success(r: FailedCallCacheDiffResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
                  case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
                  case Failure(e: CachedCallNotFoundException) => e.errorRequest(StatusCodes.NotFound)
                  case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
                }
              case Invalid(errors) =>
                val e = AggregatedMessageException("Wrong parameters for call cache diff query", errors.toList)
                e.errorRequest(StatusCodes.BadRequest)
            }
          }
        }
      }
    } ~
    path("workflows" / Segment / Segment / "timing") { (_, possibleWorkflowId) =>
      instrumentRequest {
        onComplete(validateWorkflowId(possibleWorkflowId)) {
          case Success(_) => getFromResource("workflowTimings/workflowTimings.html")
          case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
        }
      }
    } ~
    path("workflows" / Segment / Segment / "abort") { (_, possibleWorkflowId) =>
      post {
        instrumentRequest {
          def sendAbort(workflowId: WorkflowId): Route = {
            val response = workflowStoreActor.ask(WorkflowStoreActor.AbortWorkflowCommand(workflowId)).mapTo[AbortResponse]

            onComplete(response) {
              case Success(WorkflowAbortingResponse(id, restarted)) =>
                workflowManagerActor ! AbortWorkflowCommand(id, restarted)
                complete(ToResponseMarshallable(WorkflowAbortResponse(id.toString, WorkflowAborting.toString)))
              case Success(WorkflowAbortFailureResponse(_, e: IllegalStateException)) =>
                /*
                  Note that this is currently impossible to reach but was left as-is during the transition to akka http.
                  When aborts get fixed, this should be looked at.
                */
                e.errorRequest(StatusCodes.Forbidden)
              case Success(WorkflowAbortFailureResponse(_, e: WorkflowNotFoundException)) => e.errorRequest(StatusCodes.NotFound)
              case Success(WorkflowAbortFailureResponse(_, e)) => e.errorRequest(StatusCodes.InternalServerError)
              case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
              case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
            }
          }

          Try(WorkflowId.fromString(possibleWorkflowId)) match {
            case Success(workflowId) => sendAbort(workflowId)
            case Failure(_) => InvalidWorkflowException(possibleWorkflowId).failRequest(StatusCodes.BadRequest)
          }
        }
      }
    } ~
    path("workflows" / Segment / Segment / "labels") { (_, possibleWorkflowId) =>
      get {
        instrumentRequest { metadataBuilderRequest(possibleWorkflowId, (w: WorkflowId) => GetLabels(w)) }
      } ~
      patch {
        entity(as[Map[String, String]]) { parameterMap =>
          instrumentRequest {
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
      }
    } ~
  path("workflows" / Segment) { _ =>
      post {
        instrumentRequest {
          entity(as[Multipart.FormData]) { formData =>
            submitRequest(formData, isSingleSubmission = true)
          }
        }
      }
    } ~
  path("workflows" / Segment / "batch") { _ =>
    post {
      instrumentRequest {
        entity(as[Multipart.FormData]) { formData =>
          submitRequest(formData, isSingleSubmission = false)
        }
      }
    }
  }

  private def submitRequest(formData: Multipart.FormData, isSingleSubmission: Boolean): Route = {
    val allParts: Future[Map[String, ByteString]] = formData.parts.mapAsync[(String, ByteString)](1) {
      bodyPart => bodyPart.toStrict(duration).map(strict => bodyPart.name -> strict.entity.data)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)

    def toResponse(workflowId: WorkflowId): WorkflowSubmitResponse = {
      WorkflowSubmitResponse(workflowId.toString, WorkflowSubmitted.toString)
    }

    def askSubmit(command: WorkflowStoreActor.WorkflowStoreActorSubmitCommand, warnings: Seq[String]): Route = {
      // NOTE: Do not blindly copy the akka-http -to- ask-actor pattern below without knowing the pros and cons.
      onComplete(workflowStoreActor.ask(command).mapTo[WorkflowStoreSubmitActor.WorkflowStoreSubmitActorResponse]) {
        case Success(w) =>
          w match {
            case WorkflowStoreSubmitActor.WorkflowSubmittedToStore(workflowId) =>
              completeResponse(StatusCodes.Created, toResponse(workflowId), warnings)
            case WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore(workflowIds) =>
              completeResponse(StatusCodes.Created, workflowIds.toList map toResponse, warnings)
            case WorkflowStoreSubmitActor.WorkflowSubmitFailed(throwable) =>
              throwable.failRequest(StatusCodes.BadRequest, warnings)
          }
        case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
        case Failure(e) =>
          e.failRequest(StatusCodes.InternalServerError, warnings)
      }
    }

    onComplete(allParts) {
      case Success(data) =>
        PartialWorkflowSources.fromSubmitRoute(data, allowNoInputs = isSingleSubmission) match {
          case Success(workflowSourceFiles) if isSingleSubmission && workflowSourceFiles.size == 1 =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            askSubmit(WorkflowStoreActor.SubmitWorkflow(workflowSourceFiles.head), warnings)
          // Catches the case where someone has gone through the single submission endpoint w/ more than one workflow
          case Success(workflowSourceFiles) if isSingleSubmission =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            val e = new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint.")
            e.failRequest(StatusCodes.BadRequest, warnings)
          case Success(workflowSourceFiles) =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            askSubmit(
              WorkflowStoreActor.BatchSubmitWorkflows(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList)),
              warnings)
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
          case UnrecognizedWorkflowId => throw UnrecognizedWorkflowException(w)
          case FailedToCheckWorkflowId(t) => throw t
        }
      case Failure(_) => Future.failed(InvalidWorkflowException(possibleWorkflowId))
    }
  }

  private def metadataBuilderRequest(possibleWorkflowId: String, request: WorkflowId => ReadAction): Route = {
    val metadataBuilderActor = actorRefFactory.actorOf(MetadataBuilderActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName)
    val response = validateWorkflowId(possibleWorkflowId) flatMap { w => metadataBuilderActor.ask(request(w)).mapTo[MetadataBuilderActorResponse] }

    onComplete(response) {
      case Success(r: BuiltMetadataResponse) => complete(r.response)
      case Success(r: FailedMetadataResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
      case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
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
          complete(ToResponseMarshallable(w.response))
        }
      case Success(w: WorkflowQueryFailure) => w.reason.failRequest(StatusCodes.BadRequest)
      case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }
}

object CromwellApiService {
  import spray.json._

  implicit class EnhancedThrowable(val e: Throwable) extends AnyVal {
    def failRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, APIResponse.fail(e).toJson.prettyPrint, warnings)
    }
    def errorRequest(statusCode: StatusCode, warnings: Seq[String] = Vector.empty): Route = {
      completeResponse(statusCode, APIResponse.error(e).toJson.prettyPrint, warnings)
    }
  }

  def completeResponse[A](statusCode: StatusCode, value: A, warnings: Seq[String])
                         (implicit mt: ToEntityMarshaller[A]): Route = {
    val warningHeaders = warnings.toIndexedSeq map { warning =>
      /*
      Need a quoted string.
      https://stackoverflow.com/questions/7886782

      Using a poor version of ~~#!
      https://github.com/akka/akka-http/blob/v10.0.9/akka-http-core/src/main/scala/akka/http/impl/util/Rendering.scala#L206
       */
      val quotedString = "\"" + warning.replaceAll("\"","\\\\\"").replaceAll("[\\r\\n]+", " ").trim + "\""

      // https://www.w3.org/Protocols/rfc2616/rfc2616-sec14.html#sec14.46
      RawHeader("Warning", s"299 cromwell/$cromwellVersion $quotedString")
    }

    complete((statusCode, warningHeaders, value))
  }

  final case class BackendResponse(supportedBackends: List[String], defaultBackend: String)

  final case class UnrecognizedWorkflowException(id: WorkflowId) extends Exception(s"Unrecognized workflow ID: $id")

  final case class InvalidWorkflowException(possibleWorkflowId: String) extends Exception(s"Invalid workflow ID: '$possibleWorkflowId'.")

  val cromwellVersion = VersionUtil.getVersion("cromwell-engine")
  val swaggerUiVersion = VersionUtil.getVersion("swagger-ui", VersionUtil.sbtDependencyVersion("swaggerUi"))
  val backendResponse = BackendResponse(BackendConfiguration.AllBackendEntries.map(_.name).sorted, BackendConfiguration.DefaultBackendEntry.name)
  val versionResponse = JsObject(Map("cromwell" -> cromwellVersion.toJson))
  val serviceShuttingDownResponse = new Exception("Cromwell service is shutting down.").failRequest(StatusCodes.ServiceUnavailable)
}
