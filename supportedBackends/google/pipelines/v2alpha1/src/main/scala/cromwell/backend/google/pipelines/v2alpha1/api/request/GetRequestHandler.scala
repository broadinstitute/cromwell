package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.time.OffsetDateTime

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.v2alpha1.model._
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.google.pipelines.common.api.RunStatus.{Initializing, Running, Success}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.Deserialization._
import cromwell.backend.google.pipelines.v2alpha1.api.request.ErrorReporter._
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.ExecutionEvent

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

  private [request] def interpretOperationStatus(operation: Operation, pollingRequest: PAPIStatusPollRequest): RunStatus = {
    require(operation != null, "Operation must not be null.")

    try {
      if (operation.getDone) {
        val metadata = operation.getMetadata.asScala.toMap
        // Deserialize the response
        val events: List[Event] = operation.events.fallBackTo(List.empty)(pollingRequest.workflowId -> operation)
        val pipeline: Option[Pipeline] = operation.pipeline.toErrorOr.fallBack(pollingRequest.workflowId -> operation)
        val workerEvent: Option[WorkerAssignedEvent] = findEvent[WorkerAssignedEvent](events).flatMap(_(pollingRequest.workflowId -> operation))
        val executionEvents = getEventList(metadata, events).toList
        val machineType = pipeline.map(_.getResources.getVirtualMachine.getMachineType)
        // preemptible is only used if the job fails, as a heuristic to guess if the VM was preempted.
        // If we can't get the value of preempted we still need to return something, returning false will not make the failure count
        // as a preemption which seems better than saying that it was preemptible when we really don't know
        val preemptible = pipeline.exists(_.getResources.getVirtualMachine.getPreemptible.booleanValue())
        val instanceName = workerEvent.map(_.getInstance())
        val zone = workerEvent.map(_.getZone)

        // If there's an error, generate an unsuccessful status. Otherwise, we were successful!
        Option(operation.getError) match {
          case Some(error) =>
            val actions = pipeline.map(_.getActions.asScala.toList).toList.flatten
            val errorReporter = new ErrorReporter(machineType, preemptible, executionEvents, zone, instanceName, actions, operation, pollingRequest.workflowId)
            errorReporter.toUnsuccessfulRunStatus(error, events)
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

    starterEvent.toList ++ events.map(toExecutionEvent)
  }
}
