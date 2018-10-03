package cromwell.perf

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
  val cacheCopyRetries: Int = {
    val numFailures: Option[Int] = for {
      cachingObject <- callCaching
      failuresMap <- cachingObject.hitFailures
    } yield failuresMap.size

    numFailures.getOrElse(0)
  }

  /***
    * @return Time (in Duration) job spent in just CheckingCallCache state
    */
  val timeInCheckingCallCacheState: Duration = {
    val eventsRelatedToCC = executionEvents.collect {
      case event if event.description.equalsIgnoreCase("CheckingCallCache") => Duration.between(event.startTime, event.endTime)
    }

    if(eventsRelatedToCC.nonEmpty) eventsRelatedToCC.reduce(_ plus _) else Duration.ZERO
  }

  /**
    * @return Time (in Duration) job spent in Call Caching states
    */
  val timeInCallCachingState: Duration = {
    val eventsRelatedToCC = executionEvents.filter(event => callCachingEventStates.exists(state => state.equalsIgnoreCase(event.description)))

    if (eventsRelatedToCC.nonEmpty) {
      val eventsSortedByStartTime = eventsRelatedToCC.sortBy(_.startTime)
      Duration.between(eventsSortedByStartTime.head.startTime, eventsSortedByStartTime.last.endTime)
    }
    else Duration.ZERO
  }

  /**
    * @return Time (in Duration) job spent in job preparation state i.e time between job submission and running it
    *         (this doesn't consider the time it spent in call caching states
    */
  val timeInJobPreparation: Duration = {
    val durationOfEventsInPreparationState = executionEvents.collect {
      case event if jobPreparationEventStates.exists(_.equalsIgnoreCase(event.description)) => Duration.between(event.startTime, event.endTime)
    }

    if(durationOfEventsInPreparationState.nonEmpty) durationOfEventsInPreparationState.reduce(_ plus _) else Duration.ZERO
  }

  /**
    *
    * @return Time (in Duration) job spent in fetching and copying cache hit(s)
    */
  val timeForFetchingAndCopyingCacheHit: Duration = {
    val durationOfEventsInCacheCopyingState = executionEvents.collect {
      case event if cacheCopyingEventStates.exists(_.equalsIgnoreCase(event.description)) => Duration.between(event.startTime, event.endTime)
    }

    if(durationOfEventsInCacheCopyingState.nonEmpty) durationOfEventsInCacheCopyingState.reduce(_ plus _) else Duration.ZERO
  }
}
