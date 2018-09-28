package perf

import java.time.{Duration, OffsetDateTime}

import io.circe.JsonObject

case class CallCaching(hit: Option[Boolean],
                       result: Option[String],
                       hitFailures: Option[Seq[Map[String, Seq[JsonObject]]]])


case class ExecutionEvents(startTime: OffsetDateTime,
                           description: String,
                           endTime: OffsetDateTime)


case class Call(shardIndex: Int,
                start: OffsetDateTime,
                end: OffsetDateTime,
                callCaching: Option[CallCaching],
                executionEvents: Seq[ExecutionEvents]) {
  val callCachingEventStates = List("CheckingCallCache", "FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")
  val jobPreparationEventStates = List("Pending", "RequestingExecutionToken", "WaitingForValueStore", "PreparingJob", "CheckingJobStore")
  val cacheCopyingEventStates = List("FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")

  /**
    * @return Cache copy retries before getting a successful hit or running the job
    */
  def cacheCopyRetries: Int = {
    val hitFailuresMap = callCaching.map(callCachingObject => callCachingObject.hitFailures.getOrElse(Map.empty))

    if (hitFailuresMap.isDefined) hitFailuresMap.size else 0
  }

  /**
    * @return Time (in Duration) job spent in Call Caching states
    */
  def timeInCallCachingState: Duration = {
    val eventsRelatedToCC = executionEvents.filter(event => callCachingEventStates.exists(state => state.equalsIgnoreCase(event.description)))

    if (eventsRelatedToCC.nonEmpty) {
      val eventsSortedByStartTime = eventsRelatedToCC.sortBy(x => x.startTime)
      Duration.between(eventsSortedByStartTime.head.startTime, eventsSortedByStartTime.last.endTime)
    }
    else Duration.ZERO
  }

  /**
    * @return Time (in Duration) job spent in job preparation state i.e time between job submission and running it
    *         (this doesn't consider the time it spent in call caching states
    */
  def timeInJobPreparation: Duration = {
    val durationOfEventsInPreparationState = executionEvents.filter(event => jobPreparationEventStates.exists(state => state.equalsIgnoreCase(event.description)))
      .map(event => Duration.between(event.startTime, event.endTime))

    if(durationOfEventsInPreparationState.nonEmpty) durationOfEventsInPreparationState.reduce(_ plus _) else Duration.ZERO
  }

  /**
    *
    * @return Time (in Duration) job spent in fetching and copying cache hit(s)
    */
  def timeForFetchingAndCopyingCacheHit: Duration = {
    val durationOfEventsInCacheCopyingState = executionEvents.filter(event => cacheCopyingEventStates.exists(state => state.equalsIgnoreCase(event.description)))
      .map(event => Duration.between(event.startTime, event.endTime))

    if(durationOfEventsInCacheCopyingState.nonEmpty) durationOfEventsInCacheCopyingState.reduce(_ plus _) else Duration.ZERO
  }
}
