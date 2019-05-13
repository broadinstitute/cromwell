package cromwell.backend.google.pipelines.v1alpha2.api.request

import java.time.OffsetDateTime
import java.util.{ArrayList => JArrayList}

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonErrorContainer}
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import com.google.api.client.util.{ArrayMap => GArrayMap}
import com.google.api.services.genomics.model.Operation
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.google.pipelines.common.api.RunStatus.{Initializing, Running, Success, UnsuccessfulRunStatus}
import cromwell.core.ExecutionEvent
import io.grpc.Status

import scala.collection.JavaConverters._
import scala.concurrent.{Future, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Try, Success => TrySuccess}

trait GetRequestHandler { this: RequestHandler =>
  private [api] def statusPollResultHandler(originalRequest: PAPIStatusPollRequest, completionPromise: Promise[Try[Unit]], pollingManager: ActorRef) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      originalRequest.requester ! interpretOperationStatus(operation)
      completionPromise.trySuccess(TrySuccess(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      pollingManager ! PipelinesApiStatusQueryFailed(originalRequest, new SystemPAPIApiException(GoogleJsonException(e, responseHeaders)))
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  def enqueueStatusPollInBatch(pollingRequest: PAPIStatusPollRequest, batch: BatchRequest, pollingManager: ActorRef): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = statusPollResultHandler(pollingRequest, completionPromise, pollingManager)
    addStatusPollToBatch(pollingRequest.httpRequest, batch, resultHandler)
    completionPromise.future
  }

  // Done so that the specs can override this behavior
  private [api] def interpretOperationStatus(op: Operation) = GetRequestHandler.interpretOperationStatus(op)

  private def addStatusPollToBatch(request: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]): Unit = {
    /*
      * Manually enqueue the request instead of doing it through the RunPipelineRequest
      * as it would unnecessarily rebuild the request (which we already have)
    */
    batch.queue(request, classOf[Operation], classOf[GoogleJsonErrorContainer], resultHandler)
    ()
  }
}

private[api] object GetRequestHandler {

  private val AcceptableEvents = Set("start", "pulling-image", "localizing-files", "running-docker", "delocalizing-files", "ok", "fail", "start-shutdown", "preempted")

  implicit class RunOperationExtension(val operation: Operation) extends AnyVal {
    def hasStarted = operation.getMetadata.asScala.get("startTime") isDefined
  }

  def interpretOperationStatus(op: Operation): RunStatus = {
    require(op != null, "Operation must not be null.")
    try {
      if (op.getDone) {
        lazy val eventList = getEventList(op)
        lazy val computeEngineOption = for {
          runtimeMetadata <- op.getMetadata.asScala.get("runtimeMetadata")
          computeEngine <- runtimeMetadata.asInstanceOf[GArrayMap[String, Object]].asScala.get("computeEngine")
        } yield computeEngine.asInstanceOf[GArrayMap[String, String]].asScala

        lazy val machineType = computeEngineOption.flatMap(_.get("machineType"))
        lazy val instanceName = computeEngineOption.flatMap(_.get("instanceName"))
        lazy val zone = computeEngineOption.flatMap(_.get("zone"))

        val preemptible = {
          for {
            request <- op.getMetadata.asScala.get("request")
            pipelineArgs <- request.asInstanceOf[GArrayMap[String, Object]].asScala.get("pipelineArgs")
            resources <- pipelineArgs.asInstanceOf[GArrayMap[String, Object]].asScala.get("resources")
            preemptible <- resources.asInstanceOf[GArrayMap[String, Object]].asScala.get("preemptible")
          } yield preemptible.asInstanceOf[Boolean]
        } getOrElse false

        // If there's an error, generate an unsuccessful status. Otherwise, we were successful!
        Option(op.getError) match {
          case Some(error) =>
            val errorCode = Status.fromCodeValue(error.getCode)
            UnsuccessfulRunStatus(errorCode, Option(error.getMessage), eventList, machineType, zone, instanceName, preemptible)
          case None => Success(eventList, machineType, zone, instanceName)
        }
      } else if (op.hasStarted) {
        Running
      } else {
        Initializing
      }
    } catch {
      case npe: NullPointerException =>
        throw new RuntimeException(s"Caught NPE while processing operation ${op.getName}: $op", npe)
    }
  }

  def getEventList(op: Operation): Seq[ExecutionEvent] = {
    val metadata = op.getMetadata.asScala.toMap

    val starterEvents: Seq[ExecutionEvent] = Seq(
      eventIfExists("createTime", metadata, "waiting for quota"),
      eventIfExists("startTime", metadata, "initializing VM")).flatten

    val eventsList = for {
      events <- metadata.get("events").toSeq
      entry <- events.asInstanceOf[JArrayList[GArrayMap[String, String]]].asScala
    } yield ExecutionEvent(entry.get("description"), OffsetDateTime.parse(entry.get("startTime")))

    val filteredEventsList: Seq[ExecutionEvent] = eventsList filter { i => AcceptableEvents.contains(i.name) }

    // A little bit ugly... the endTime of the jes operation can actually be before the final "event" time, due to differences
    // in the reported precision. As a result, we have to make sure it all lines up nicely:
    val finalEvent = getCromwellPollIntervalEvent(metadata, filteredEventsList)

    starterEvents ++ filteredEventsList :+ finalEvent
  }

  private def getCromwellPollIntervalEvent(metadata: Map[String, AnyRef], eventsList: Seq[ExecutionEvent]) = {
    {
      val jesReportedEndTime = eventIfExists("endTime", metadata, "cromwell poll interval")
      val finalEventsListTime = if (eventsList.nonEmpty) Some(eventsList.last.offsetDateTime) else None

      (jesReportedEndTime, finalEventsListTime) match {
        case (Some(jesEndTime), Some(finalEventTime)) =>
          if (jesEndTime.offsetDateTime isAfter finalEventTime) jesEndTime else jesEndTime.copy(offsetDateTime = finalEventTime)
        case (Some(jesEndTime), None) => jesEndTime
        case (None, Some(finalEventTime)) => ExecutionEvent("cromwell poll interval", finalEventTime)
        case (None, None) =>
          throw new IllegalArgumentException("Both jesReportedEndTime and finalEventsListTime were None.")
      }
    }
  }

  private def eventIfExists(name: String, metadata: Map[String, AnyRef], eventName: String): Option[ExecutionEvent] = {
    metadata.get(name) map { time => ExecutionEvent(eventName, OffsetDateTime.parse(time.toString)) }
  }
}
