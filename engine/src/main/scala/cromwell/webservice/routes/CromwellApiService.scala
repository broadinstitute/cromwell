package cromwell.webservice.routes

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
import cromwell.core.abort.{AbortResponse, WorkflowAbortFailureResponse, WorkflowAbortedResponse, WorkflowAbortingResponse}
import cromwell.core.{path => _, _}
import cromwell.engine.backend.BackendConfiguration
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.engine.workflow.WorkflowManagerActor.{AbortWorkflowCommand, WorkflowNotFoundException}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffActor.{BuiltCallCacheDiffResponse, CachedCallNotFoundException, CallCacheDiffActorResponse, FailedCallCacheDiffResponse}
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheDiffActor, CallCacheDiffQueryParameter}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.NotInOnHoldStateException
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreEngineActor, WorkflowStoreSubmitActor}
import cromwell.server.CromwellShutdown
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{GetCurrentStatus, StatusCheckResponse}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, TimeoutException}
import scala.util.{Failure, Success, Try}

trait CromwellApiService extends HttpInstrumentation with MetadataRouteSupport {
  import CromwellApiService._

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
        completeResponse(StatusCodes.Forbidden, APIResponse.fail(new RuntimeException("The /stats endpoint is currently disabled.")), warnings = Seq.empty)
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
        case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
        case Failure(_) => new RuntimeException("Unable to gather engine status").failRequest(StatusCodes.InternalServerError)
      }
    }
  )

  val workflowRoutes =
    path("workflows" / Segment / "backends") { _ =>
      get { instrumentRequest { complete(ToResponseMarshallable(backendResponse)) } }
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
                  case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
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
        onComplete(validateWorkflowId(possibleWorkflowId, serviceRegistryActor)) {
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
              case Success(WorkflowAbortedResponse(id)) =>
                complete(ToResponseMarshallable(WorkflowAbortResponse(id.toString, WorkflowAborted.toString)))
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
              case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
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
  } ~
  path("workflows" / Segment / Segment / "releaseHold") { (_, possibleWorkflowId) =>
    post {
      instrumentRequest {
        val response = validateWorkflowId(possibleWorkflowId, serviceRegistryActor) flatMap { workflowId =>
          workflowStoreActor.ask(WorkflowStoreActor.WorkflowOnHoldToSubmittedCommand(workflowId)).mapTo[WorkflowStoreEngineActor.WorkflowOnHoldToSubmittedResponse]
        }
        onComplete(response){
          case Success(WorkflowStoreEngineActor.WorkflowOnHoldToSubmittedFailure(_, e: NotInOnHoldStateException)) => e.errorRequest(StatusCodes.Forbidden)
          case Success(WorkflowStoreEngineActor.WorkflowOnHoldToSubmittedFailure(_, e)) => e.errorRequest(StatusCodes.InternalServerError)
          case Success(r: WorkflowStoreEngineActor.WorkflowOnHoldToSubmittedSuccess) => completeResponse(StatusCodes.OK, toResponse(r.workflowId, WorkflowSubmitted), Seq.empty)
          case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
          case Failure(e: InvalidWorkflowException) => e.failRequest(StatusCodes.BadRequest)
          case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
        }
      }
    }
  } ~ metadataRoutes

  private def toResponse(workflowId: WorkflowId, workflowState: WorkflowState): WorkflowSubmitResponse = {
    WorkflowSubmitResponse(workflowId.toString, workflowState.toString)
  }

  private def submitRequest(formData: Multipart.FormData, isSingleSubmission: Boolean): Route = {
    val allParts: Future[Map[String, ByteString]] = formData.parts.mapAsync[(String, ByteString)](1) {
      bodyPart => bodyPart.toStrict(duration).map(strict => bodyPart.name -> strict.entity.data)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)

    def getWorkflowState(workflowOnHold: Boolean): WorkflowState = {
      if (workflowOnHold)
        WorkflowOnHold
      else WorkflowSubmitted
    }

    def askSubmit(command: WorkflowStoreActor.WorkflowStoreActorSubmitCommand, warnings: Seq[String], workflowState: WorkflowState): Route = {
      // NOTE: Do not blindly copy the akka-http -to- ask-actor pattern below without knowing the pros and cons.
      onComplete(workflowStoreActor.ask(command).mapTo[WorkflowStoreSubmitActor.WorkflowStoreSubmitActorResponse]) {
        case Success(w) =>
          w match {
            case WorkflowStoreSubmitActor.WorkflowSubmittedToStore(workflowId, _) =>
              completeResponse(StatusCodes.Created, toResponse(workflowId, workflowState), warnings)
            case WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore(workflowIds, _) =>
              completeResponse(StatusCodes.Created, workflowIds.toList.map(toResponse(_, workflowState)), warnings)
            case WorkflowStoreSubmitActor.WorkflowSubmitFailed(throwable) =>
              throwable.failRequest(StatusCodes.BadRequest, warnings)
          }
        case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
        case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
        case Failure(e) => e.failRequest(StatusCodes.InternalServerError, warnings)
      }
    }

    onComplete(allParts) {
      case Success(data) =>
        PartialWorkflowSources.fromSubmitRoute(data, allowNoInputs = isSingleSubmission) match {
          case Success(workflowSourceFiles) if isSingleSubmission && workflowSourceFiles.size == 1 =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            askSubmit(WorkflowStoreActor.SubmitWorkflow(workflowSourceFiles.head), warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          // Catches the case where someone has gone through the single submission endpoint w/ more than one workflow
          case Success(workflowSourceFiles) if isSingleSubmission =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            val e = new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint.")
            e.failRequest(StatusCodes.BadRequest, warnings)
          case Success(workflowSourceFiles) =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            askSubmit(
              WorkflowStoreActor.BatchSubmitWorkflows(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList)),
              warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          case Failure(t) => t.failRequest(StatusCodes.BadRequest)
        }
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
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

  def validateWorkflowId(possibleWorkflowId: String,
                         serviceRegistryActor: ActorRef)
                        (implicit timeout: Timeout, executor: ExecutionContext): Future[WorkflowId] = {
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
