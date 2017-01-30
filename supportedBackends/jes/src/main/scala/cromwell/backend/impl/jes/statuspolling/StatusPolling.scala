package cromwell.backend.impl.jes.statuspolling

import java.time.OffsetDateTime
import java.util.{ArrayList => JArrayList}

import com.google.api.client.util.{ArrayMap => GArrayMap}
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.RunStatus._
import cromwell.backend.impl.jes.{JesAsyncBackendJobExecutionActor, Run, RunStatus}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{GoogleJsonException, JesApiException, JesApiQueryFailed, JesStatusPollQuery}
import cromwell.core.ExecutionEvent

import scala.language.postfixOps
import scala.collection.JavaConverters._
import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Try, Success => TrySuccess}

private[statuspolling] trait StatusPolling { this: JesPollingActor =>

  private def statusPollResultHandler(originalRequest: JesStatusPollQuery, completionPromise: Promise[Try[Unit]]) = new JsonBatchCallback[Operation] {
    override def onSuccess(operation: Operation, responseHeaders: HttpHeaders): Unit = {
      originalRequest.requester ! interpretOperationStatus(operation)
      completionPromise.trySuccess(TrySuccess(()))
      ()
    }

    override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
      pollingManager ! JesApiQueryFailed(originalRequest, new JesApiException(GoogleJsonException(e, responseHeaders)))
      completionPromise.trySuccess(Failure(new Exception(mkErrorString(e))))
      ()
    }
  }

  def enqueueStatusPollInBatch(pollingRequest: JesStatusPollQuery, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = statusPollResultHandler(pollingRequest, completionPromise)
    addStatusPollToBatch(pollingRequest.run, batch, resultHandler)
    completionPromise.future
  }

  // Done so that the specs can override this behaviour
  def interpretOperationStatus(op: Operation) = StatusPolling.interpretOperationStatus(op)

  def addStatusPollToBatch(run: Run, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]) =
    run.getOperationCommand.queue(batch, resultHandler)
}

private[statuspolling] object StatusPolling {

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

        // If there's an error, generate a Failed status. Otherwise, we were successful!
        Option(op.getError) match {
          case None => Success(eventList, machineType, zone, instanceName)
          case Some(error) if error.getCode == JesAsyncBackendJobExecutionActor.JesPreemption =>
            Preempted(error.getCode, Option(error.getMessage), eventList, machineType, zone, instanceName)
          case Some(error) =>
            Failed(error.getCode, Option(error.getMessage), eventList, machineType, zone, instanceName)
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
