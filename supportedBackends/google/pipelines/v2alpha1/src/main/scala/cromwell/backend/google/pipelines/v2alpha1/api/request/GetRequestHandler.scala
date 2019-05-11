package cromwell.backend.google.pipelines.v2alpha1.api.request

import java.time.OffsetDateTime

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.v2alpha1.model._
import common.validation.Validation._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.google.pipelines.common.api.RunStatus.{Initializing, Running, Success, UnsuccessfulRunStatus}
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.Key
import cromwell.backend.google.pipelines.v2alpha1.api.Deserialization._
import cromwell.backend.google.pipelines.v2alpha1.api.request.ErrorReporter._
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.ExecutionEvent
import io.grpc.Status
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
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
        case TrySuccess(googleError) => new SystemPAPIApiException(GoogleJsonException(googleError, response.getHeaders))
        case Failure(_) => new SystemPAPIApiException(new RuntimeException(s"Failed to get status for operation ${pollingRequest.jobId.jobId}: HTTP Status Code: ${response.getStatusCode}"))
      }
      pollingManager ! PipelinesApiStatusQueryFailed(pollingRequest, failure)
      Failure(failure)
  } recover {
    case e =>
      pollingManager ! PipelinesApiStatusQueryFailed(pollingRequest, new SystemPAPIApiException(e))
      Failure(e)
  }

  private [request] def interpretOperationStatus(operation: Operation, pollingRequest: PAPIStatusPollRequest): RunStatus = {
    if (Option(operation).isEmpty) {
      // It is possible to receive a null via an HTTP 200 with no response. If that happens, handle it and don't crash.
      // https://github.com/googleapis/google-http-java-client/blob/v1.28.0/google-http-client/src/main/java/com/google/api/client/http/HttpResponse.java#L456-L458
      val errorMessage = "Operation returned as empty"
      UnsuccessfulRunStatus(Status.UNKNOWN, Option(errorMessage), Nil, None, None, None, wasPreemptible = false)
    } else {
      try {
        if (operation.getDone) {
          val metadata = operation.metadata
          // Deserialize the response
          val events: List[Event] = operation.events.fallBackTo(List.empty)(pollingRequest.workflowId -> operation)
          val pipeline: Option[Pipeline] = operation.pipeline.flatMap(
            _.toErrorOr.fallBack(pollingRequest.workflowId -> operation)
          )
          val actions: List[Action] = pipeline
            .flatMap(pipelineValue => Option(pipelineValue.getActions))
            .map(_.asScala)
            .toList
            .flatten
          val workerEvent: Option[WorkerAssignedEvent] =
            findEvent[WorkerAssignedEvent](events).flatMap(_ (pollingRequest.workflowId -> operation))
          val executionEvents = getEventList(metadata, events, actions)
          val virtualMachineOption = for {
            pipelineValue <- pipeline
            resources <- Option(pipelineValue.getResources)
            virtualMachine <- Option(resources.getVirtualMachine)
          } yield virtualMachine
          // Correlate `executionEvents` to `actions` to potentially assign a grouping into the appropriate events.
          val machineType = virtualMachineOption.flatMap(virtualMachine => Option(virtualMachine.getMachineType))
          /*
          preemptible is only used if the job fails, as a heuristic to guess if the VM was preempted.
          If we can't get the value of preempted we still need to return something, returning false will not make the
          failure count as a preemption which seems better than saying that it was preemptible when we really don't know
           */
          val preemptibleOption = for {
            pipelineValue <- pipeline
            resources <- Option(pipelineValue.getResources)
            virtualMachine <- Option(resources.getVirtualMachine)
            preemptible <- Option(virtualMachine.getPreemptible)
          } yield preemptible
          val preemptible = preemptibleOption.exists(_.booleanValue)
          val instanceName = workerEvent.flatMap(workerAssignedEvent => Option(workerAssignedEvent.getInstance()))
          val zone = workerEvent.flatMap(workerAssignedEvent => Option(workerAssignedEvent.getZone))

          // If there's an error, generate an unsuccessful status. Otherwise, we were successful!
          Option(operation.getError) match {
            case Some(error) =>
              val errorReporter = new ErrorReporter(
                machineType,
                preemptible,
                executionEvents,
                zone,
                instanceName,
                actions,
                operation,
                pollingRequest.workflowId
              )
              errorReporter.toUnsuccessfulRunStatus(error, events)
            case None => Success(executionEvents, machineType, zone, instanceName)
          }
        } else if (operation.hasStarted) {
          Running
        } else {
          Initializing
        }
      } catch {
        case nullPointerException: NullPointerException =>
          throw new RuntimeException(
            s"Caught NPE while interpreting operation ${operation.getName}: " +
              s"${ExceptionUtils.getStackTrace(nullPointerException)}. " +
              s"JSON was $operation",
            nullPointerException
          )
      }
    }
  }

  private def getEventList(metadata: Map[String, AnyRef], events: List[Event], actions: List[Action]): List[ExecutionEvent] = {
    val starterEvent: Option[ExecutionEvent] = {
      metadata.get("createTime") map { time => ExecutionEvent("waiting for quota", OffsetDateTime.parse(time.toString)) }
    }

    // Map action indexes to event types. Action indexes are 1-based for some reason.
    val actionIndexToEventType: Map[Int, String] = List(Key.Logging, Key.Tag).flatMap { k =>
      actions.zipWithIndex collect { case (a, i) if a.getLabels.containsKey(k) => (i + 1) -> a.getLabels.get(k) } } toMap

    val executionEvents = events.map(toExecutionEvent(actionIndexToEventType))
    // The Docker image used for CWL output parsing causes some complications for the timing diagram. Docker image
    // pulling is done automatically by PAPI v2 and therefore does not correspond to any Actions generated in Cromwell.
    // Since there are no Actions there are no labeled Actions and therefore Docker pull events do not get grouped into
    // any labeled categories. For regular prefetched images this is okay since they don't conflict with any other labeled
    // activity and we actually want to see the details of those pulls on the timing diagram. But the CWL output parsing
    // image pull happens in the middle of other 'Delocalization' events. Since it is not a 'Delocalization' event it appears
    // to the timing diagram logic to be concurrent to 'Delocalization' and another row is wrongly added to the timing diagram.
    // The logic here checks if the event description matches a CWL Docker image pull and if so suppresses publication,
    // which is enough to make the timing diagram happy.
    val startDelocalization: Option[ExecutionEvent] =
      executionEvents filter { _.grouping.contains("Delocalization") } sortBy (_.offsetDateTime) headOption

    val filteredExecutionEvents = startDelocalization match {
      case None => executionEvents // Can't do filtering without a start time for Delocalization.
      case Some(start) =>
        executionEvents filterNot { e => (e.name.startsWith("Started pulling ") || e.name.startsWith("Stopped pulling ")) && e.offsetDateTime.compareTo(start.offsetDateTime) > 0 }
    }

    starterEvent.toList ++ filteredExecutionEvents
  }
}
