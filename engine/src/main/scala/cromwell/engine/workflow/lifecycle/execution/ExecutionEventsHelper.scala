package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import cromwell.core.{ExecutionEvent, JobKey}
import cromwell.engine.workflow.lifecycle.execution.ExecutionEventsHelper.ExecutionEventMark
import cromwell.services.metadata.CallMetadataKeys
import cromwell.services.metadata.MetadataService.PutMetadataAction

object ExecutionEventsHelper {
  private case class ExecutionEventMark(executionEvent: ExecutionEvent, key: String)
}

trait ExecutionEventsHelper { this: CallMetadataHelper =>
  /*
   * Mark the last execution event for which we didn't have an end time (and set it to null).
   * Next event that comes in we can override the mark event's end time with the new event's start time.
   */
  private var lastExecutionEventMark: Option[ExecutionEventMark] = None

  /**
    * Push several execution events. End times will be all known except for the last event.
    */
  def pushExecutionEventsToMetadataService(jobKey: JobKey, events: Seq[ExecutionEvent]) = {
    val sorted = events.sortBy(_.offsetDateTime)
    sorted.sliding(2).foreach({
      case Seq(event1, event2) => pushExecutionEventToMetadataService(jobKey)(event1, Option(event2.offsetDateTime))
    })

    // Push the last event with an unknown end time
    sorted.lastOption.foreach(pushExecutionEventToMetadataService(jobKey)(_))
  }

  /**
    * Push an execution event with an optional end time.
    */
  def pushExecutionEventToMetadataService(jobKey: JobKey)(event: ExecutionEvent, endTime: Option[OffsetDateTime] = None) = {
    val eventKey = s"${CallMetadataKeys.ExecutionEvents}[$randomNumberString]"

    // Override the last endTime if we have a mark, and its start time is before this event
    val overrideLastEndDate = lastExecutionEventMark match {
      case Some(ExecutionEventMark(markEvent, markEventKey)) if markEvent.offsetDateTime.isBefore(event.offsetDateTime) => 
        List(metadataEvent(jobKey, s"$markEventKey:endTime", event.offsetDateTime))
      case _ => List.empty
    }
    
    val events = overrideLastEndDate ++ List(
      metadataEvent(jobKey, s"$eventKey:description", event.name),
      metadataEvent(jobKey, s"$eventKey:startTime", event.offsetDateTime),
      metadataEvent(jobKey, s"$eventKey:endTime", endTime.orNull)
    )

    // If the end time is unknown, record the key so we can override its end time on the next event
    lastExecutionEventMark = if (endTime.isEmpty) Option(ExecutionEventMark(event, eventKey)) else None

    serviceRegistryActor ! PutMetadataAction(events)
  }

  /**
    * If an event was left with an unknown end time, set the end time to now.
    */
  def terminateExecutionEvents(jobKey: JobKey) = {
    val events = lastExecutionEventMark.map(marker => metadataEvent(jobKey, s"${marker.key}:endTime", OffsetDateTime.now())).toList
    lastExecutionEventMark = None
    serviceRegistryActor ! PutMetadataAction(events)
  }
}
