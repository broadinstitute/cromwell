package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cats.data.EitherT
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.traverse._
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.model.UnexpectedExitStatusEvent
import com.google.api.services.genomics.v2alpha1.model._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.google.pipelines.common.api.RunStatus.{Initializing, Running, Success, UnsuccessfulRunStatus}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.Deserialization._
import cromwell.backend.google.pipelines.v2alpha1.api.request.RequestHandler._
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.ExecutionEvent
import io.grpc.Status
import mouse.all._

import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Try, Success => TrySuccess}

trait GetRequestHandler { this: RequestHandler =>
  // the Genomics batch endpoint doesn't seem to be able to handle get requests on V2 operations at the moment
  // For now, don't batch the request and execute it on its own 
  def handleRequest(pollingRequest: PAPIStatusPollRequest, batch: BatchRequest, pollingManager: ActorRef)(implicit ec: ExecutionContext): Future[Try[Unit]] = Future(pollingRequest.httpRequest.execute()) map {
    case response if response.isSuccessStatusCode =>
      val operation = response.parseAs(classOf[Operation])
      pollingRequest.requester ! interpretOperationStatus(operation, pollingRequest)
      TrySuccess(())
    case response =>
      val failure = Try(GoogleJsonError.parse(GoogleAuthMode.jsonFactory, response)) match {
        case TrySuccess(googleError) => new PAPIApiException(GoogleJsonException(googleError, response.getHeaders))
        case Failure(_) => new PAPIApiException(new RuntimeException(s"Failed to get status for operation ${pollingRequest.jobId.jobId}: HTTP Status Code: ${response.getStatusCode}"))
      }
      pollingManager ! PipelinesApiStatusQueryFailed(pollingRequest, failure)
      Failure(failure)
  } recover {
    case e =>
      pollingManager ! PipelinesApiStatusQueryFailed(pollingRequest, new PAPIApiException(e))
      Failure(e)
  }

  private def findWorkerEvent(events: List[Event]) = EitherT.apply(
    events
      .map(_.details[WorkerAssignedEvent])
      .collectFirst({ case Some(defined) => defined.toChecked })
  )

  private [request] def interpretOperationStatus(operation: Operation, pollingRequest: PAPIStatusPollRequest): RunStatus = {
    require(operation != null, "Operation must not be null.")

    def withDefault[T](default: T, context: String)(value: ErrorOr[T]): T = value.valueOr({ errors =>
      val errorString = errors.toList.mkString(", ")
      logger.error(s"Could not $context from the pipelines API response for workflow ${pollingRequest.workflowId}, operation ${operation.getName}: $errorString")
      default
    })

    try {
      if (operation.getDone) {
        val metadata = operation.getMetadata.asScala.toMap
        // Deserialize the response
        val events: ErrorOr[List[Event]] = operation.events
        val pipeline: ErrorOr[Pipeline] = operation.pipeline.toErrorOr

        val workerEvent: Option[WorkerAssignedEvent] =
          (for {
            validEvents <- EitherT.fromEither[Option](events.toEither)
            maybeWorkerEvent <- findWorkerEvent(validEvents)
          } yield maybeWorkerEvent)
            .value.map(_.toValidated)
            .sequence[ErrorOr, WorkerAssignedEvent] |> withDefault(None, "parse details of the worker event")

        val executionEvents = events.map(getEventList(metadata, _)) |> withDefault(List.empty, "parse the events")

        // preemptible is only used if the job fails, as a heuristic to guess if the VM was preempted.
        // If we can't get the value of preempted we still need to return something, returning false will not make the failure count
        // as a preemption which seems better than saying that it was preemptible when we really don't know
        val (machineType, preemptible): (Option[String], Boolean) = {
          pipeline.map({ pipeline =>
            val machineType = pipeline.getResources.getVirtualMachine.getMachineType
            val preemptible = pipeline.getResources.getVirtualMachine.getPreemptible
            Option(machineType) -> preemptible.booleanValue()
          }) |> withDefault[(Option[String], Boolean)]((None, false), "parse the pipeline")
        }

        val instanceName = workerEvent.map(_.getInstance())
        val zone = workerEvent.map(_.getZone)

        // If there's an error, generate an unsuccessful status. Otherwise, we were successful!
        Option(operation.getError) match {
          case Some(error) =>
            val errorCode = Status.fromCodeValue(error.getCode)
            val actions: ErrorOr[List[Action]] = pipeline.map(_.getActions.asScala.toList)

            val errorMessage = (events, actions) mapN {
              case (e, a) => augmentedErrorMessage(e, a, error.getMessage)
            } getOrElse error.getMessage

            UnsuccessfulRunStatus(errorCode, Option(errorMessage), executionEvents, machineType, zone, instanceName, preemptible)
          case None => Success(executionEvents, machineType, zone, instanceName)
        }
      } else if (operation.hasStarted) {
        Running
      } else {
        Initializing
      }
    } catch {
      case npe: NullPointerException =>
        throw new RuntimeException(s"Caught NPE while processing operation ${operation.getName}: $operation", npe)
    }
  }

  private def getEventList(metadata: Map[String, AnyRef], events: List[Event]): Seq[ExecutionEvent] = {
    val starterEvent: Option[ExecutionEvent] = {
      metadata.get("createTime") map { time => ExecutionEvent("waiting for quota", OffsetDateTime.parse(time.toString)) }
    }

    starterEvent.toList ++ events.map(_.toExecutionEvent)
  }

  private def augmentedErrorMessage(events: List[Event], actions: List[Action], error: String): String = {
    import cats.instances.list._
    import cats.syntax.traverse._
    import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._

    error + events
      .flatMap(_.details[UnexpectedExitStatusEvent].map(_.toErrorOr))
      .sequence[ErrorOr, UnexpectedExitStatusEvent]
      .map(_.flatMap(_.toErrorMessage(actions)).mkString(start = "\n", sep = "\n", end = ""))
      .getOrElse("")
  }
}
