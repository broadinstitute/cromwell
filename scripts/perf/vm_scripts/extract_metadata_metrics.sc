#!/usr/bin/env amm

@
import java.time.temporal.ChronoUnit

import spray.json._
import spray.json.DefaultJsonProtocol
import better.files.File
import java.time.{Duration, OffsetDateTime}

import org.apache.commons.lang3.time.DurationFormatUtils
import com.github.tototoshi.csv._

import scala.collection.immutable.ListMap

case class CallCaching(hit: Option[Boolean],
                       result: Option[String],
                       hitFailures: Option[Seq[Map[String, Seq[JsObject]]]])

case class ExecutionEvents(startTime: OffsetDateTime,
                           description: String,
                           endTime: OffsetDateTime)

case class Call(shardIndex: Int,
                start: OffsetDateTime,
                end: OffsetDateTime,
                callCaching: Option[CallCaching],
                executionEvents: Seq[ExecutionEvents])

case class Metadata(id: String,
                    workflowName: String,
                    submission: OffsetDateTime,
                    start: OffsetDateTime,
                    end: OffsetDateTime,
                    status: String,
                    calls: Option[Map[String, Seq[Call]]])

case class WorkflowMetrics(workflowId: String,
                           workflowName: String,
                           workflowStartedAfter: String,
                           workflowRunningTime: String,
                           totalJobsPerRootWf: Int,
                           avgCacheCopyRetries: Double,
                           avgTimeInCallCachingState: String,
                           avgTimeFromSubmissionToRunning: String,
                           avgTimeForFetchingCopyingCacheHit: String)

case class ScatterWidthMetrics(taskName: String,
                               scatterWidth: Int)

case class JobMetrics(shard: Int,
                      cacheCopyRetries: Int,
                      timeInCallCachingState: String,
                      timeFromSubmissionToRunning: String,
                      timeForFetchingCopyingCacheHit: String)

case class TaskMetrics(taskName: String,
                       jobMetrics: Seq[JobMetrics])


object OffsetDateTimeJsonFormatter extends DefaultJsonProtocol {
  object OffsetDateTimeFormat extends RootJsonFormat[OffsetDateTime] {
    def write(odt: OffsetDateTime) = new JsString(odt.toString)
    def read(value: JsValue) = value match {
      case JsString(string) => OffsetDateTime.parse(string)
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into an OffsetDateTime")
    }
  }
}

object MetadataJsonProtocol extends DefaultJsonProtocol {
  implicit val OffsetDateTimeJsonFormat = OffsetDateTimeJsonFormatter.OffsetDateTimeFormat

  //  implicit val callCachingHitFailuresMsgFormat: JsonFormat[CallCachingHitFailuresMsg] = lazyFormat(jsonFormat(CallCachingHitFailuresMsg, "causedBy", "message"))
  implicit val callCachingFormat = jsonFormat3(CallCaching)
  implicit val executionEventsFormat = jsonFormat3(ExecutionEvents)
  implicit val callFormat = jsonFormat5(Call)
  implicit val metadataFormat = jsonFormat7(Metadata)
}

@main
def extractMetricsFromMetadata(filePath: String): Unit = {
  import MetadataJsonProtocol._

  val callCachingEventStates = List("CheckingCallCache", "FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")
  val jobPreparationEventStates = List("Pending", "RequestingExecutionToken", "WaitingForValueStore", "PreparingJob", "CheckingJobStore")
  val cacheCopyingEventStates = List("FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")

  var totalJobsByRootWf: Int = 0
  var totalCacheTries: Double = 0
  var totalTimeInCallCaching: Duration = Duration.ZERO
  var totalTimeInJobPreparation: Duration = Duration.ZERO
  var totalTimeInCacheCopying: Duration = Duration.ZERO

  //    val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/small_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/hello_world_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/subworkflow_hello_world_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/subworkflow_with_scatter_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/cc_after_metadata.json")
  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/cc_before_metadata.json")
  //    val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/large_scatter_with_multiple_calls_metadata.json")


  val metadataFileContent = metadataFile1.contentAsString
  val workflowMetadata = metadataFileContent.parseJson.convertTo[Metadata]

  val workflowStartedAfterTime = Duration.between(workflowMetadata.submission, workflowMetadata.start)
  val workflowRunningTimeDuration = Duration.between(workflowMetadata.start, workflowMetadata.end)

  //scatter width at expansion time
  val scatterWidthMetricsOption = workflowMetadata.calls.map(callMap => callMap.map(task => {
    if(task._2.head.shardIndex != -1) ScatterWidthMetrics(task._1, task._2.size)
    else ScatterWidthMetrics(task._1, -1)
  }))

  //metrics per task, which contains job level metrics as well
  val taskLevelMetricsOption: Option[Iterable[TaskMetrics]] = workflowMetadata.calls.map(taskMap => taskMap.map(task => {
    totalJobsByRootWf += task._2.size

    val allJobMetrics: Seq[JobMetrics] = task._2.map(call => {
      //cache copy tries before finding a successful hit
      val hitFailuresMap = call.callCaching.map(callCachingObject => callCachingObject.hitFailures.getOrElse(Map.empty))

      val cacheRetries = if(hitFailuresMap.isDefined) hitFailuresMap.get.size else 0
      totalCacheTries += cacheRetries

      //time job spent in call caching states
      val eventsRelatedToCC = call.executionEvents.filter(event => callCachingEventStates.exists(state => state.equalsIgnoreCase(event.description)))

      val timeInCallCachingState = if (eventsRelatedToCC.nonEmpty) {
        val eventsSortedByStartTime = eventsRelatedToCC.sortBy(x => x.startTime)
        val durationOfCC = Duration.between(eventsSortedByStartTime.head.startTime, eventsSortedByStartTime.last.endTime)

        totalTimeInCallCaching = totalTimeInCallCaching.plus(durationOfCC)

        durationOfCC
      }
      else Duration.ZERO


      //time job spent in job preparation state i.e time between job submission and running it
      //(this doesn't consider the time it spent in call caching states
      val durationOfEventsInPreparationState = call.executionEvents
        .filter(event => jobPreparationEventStates.exists(state => state.equalsIgnoreCase(event.description)))
        .map(event => Duration.between(event.startTime, event.endTime))

      val timeInPreparationState = if(durationOfEventsInPreparationState.nonEmpty) {
        val totalTimeInPreparation = durationOfEventsInPreparationState.reduce(_ plus _)

        totalTimeInJobPreparation = totalTimeInJobPreparation.plus(totalTimeInPreparation)

        totalTimeInPreparation
      }
      else Duration.ZERO

      //time job spent in fetching and copying cache hit(s)
      val durationOfEventsInCacheCopyingState = call.executionEvents
        .filter(event => cacheCopyingEventStates.exists(state => state.equalsIgnoreCase(event.description)))
        .map(event => Duration.between(event.startTime, event.endTime))

      val timeInCopyingCache = if(durationOfEventsInCacheCopyingState.nonEmpty){
        val totalTimeCopyingCacheHits = durationOfEventsInCacheCopyingState.reduce(_ plus _)

        totalTimeInCacheCopying = totalTimeInCacheCopying.plus(totalTimeCopyingCacheHits)

        totalTimeCopyingCacheHits
      }
      else Duration.ZERO


      JobMetrics(
        shard = call.shardIndex,
        cacheCopyRetries = cacheRetries,
        timeInCallCachingState = timeInCallCachingState,
        timeInJobPreparation = timeInPreparationState,
        timeForFetchingCopyingCacheHit = timeInCopyingCache
      )
    })

    TaskMetrics(
      taskName = task._1,
      jobMetrics = allJobMetrics
    )
  }))


  println(s"Total job ct $totalJobsByRootWf")
  println(s"Total cache retries $totalCacheTries")
  println(s"Total time in CC state $totalTimeInCallCaching avg: ${totalTimeInCallCaching.dividedBy(totalJobsByRootWf)}")
  println(s"Total time in job preparation state $totalTimeInJobPreparation avg: ${totalTimeInJobPreparation.dividedBy(totalJobsByRootWf)}")
  println(s"Total time in copying cache $totalTimeInCacheCopying avg: ${totalTimeInCacheCopying.dividedBy(totalJobsByRootWf)}")


  //workflow level metrics
  val workflowMetricsObject = if(totalJobsByRootWf > 0) {
    val avgCacheRetries = totalCacheTries/totalJobsByRootWf
    val avgTimeInCallCaching = totalTimeInCallCaching.dividedBy(totalJobsByRootWf)
    val avgTimeInJobPreparationState = totalTimeInJobPreparation.dividedBy(totalJobsByRootWf)
    val avgTimeInCacheCopying = totalTimeInCacheCopying.dividedBy(totalJobsByRootWf)

    //calculate standard deviation
    //    val stdDevOfCacheRetries = taskLevelMetricsOption.map(taskLevelMetricsList =>
    //      math.sqrt ((0.0 /: taskLevelMetricsList.map(taskMetrics => {
    //        (acc, job) => acc + math.pow(job.cacheCopyRetries - avgCacheRetries, 2)
    //      }))/totalJobsByRootWf)

    val stdDevOfCacheRetries = taskLevelMetricsOption.map(taskLevelMetricsList => math.sqrt(taskLevelMetricsList.map(taskMetrics => {
      (0.0 /: taskMetrics.jobMetrics) {
        (acc, job) => acc + math.pow(job.cacheCopyRetries - avgCacheRetries, 2)
      }
    }).sum / totalJobsByRootWf))

    val stdDevOfTimeInCallCaching = taskLevelMetricsOption.map(taskLevelMetricsList => Duration.of(math.sqrt(taskLevelMetricsList.map(taskMetrics => {
      (0.0 /: taskMetrics.jobMetrics) {
        (acc, job) => acc + math.pow(job.timeInCallCachingState.minus(avgTimeInCallCaching).toMillis, 2)
      }
    }).sum / totalJobsByRootWf).toLong, ChronoUnit.MILLIS))

    val stdDevOfTimeInJobPreparation = taskLevelMetricsOption.map(taskLevelMetricsList => Duration.of(math.sqrt(taskLevelMetricsList.map(taskMetrics => {
      (0.0 /: taskMetrics.jobMetrics) {
        (acc, job) => acc + math.pow(job.timeInJobPreparation.minus(avgTimeInJobPreparationState).toMillis, 2)
      }
    }).sum / totalJobsByRootWf).toLong, ChronoUnit.MILLIS))

    val stdDevOfTimeInCacheCopying = taskLevelMetricsOption.map(taskLevelMetricsList => Duration.of(math.sqrt(taskLevelMetricsList.map(taskMetrics => {
      (0.0 /: taskMetrics.jobMetrics) {
        (acc, job) => acc + math.pow(job.timeForFetchingCopyingCacheHit.minus(avgTimeInCacheCopying).toMillis, 2)
      }
    }).sum / totalJobsByRootWf).toLong, ChronoUnit.MILLIS))

    println(s"Standard Deviation of Cache Retries: $stdDevOfCacheRetries")
    println(s"Standard Deviation of time in Call Caching: $stdDevOfTimeInCallCaching")
    println(s"Standard Deviation of time in Job Preparation: $stdDevOfTimeInJobPreparation")
    println(s"Standard Deviation of time in Cache Copying: $stdDevOfTimeInCacheCopying")


    WorkflowMetrics(
      workflowId = workflowMetadata.id,
      workflowName = workflowMetadata.workflowName,
      workflowStartedAfter = workflowStartedAfterTime,
      workflowRunningTime = workflowRunningTimeDuration,
      totalJobsPerRootWf = totalJobsByRootWf,
      avgCacheCopyRetries = avgCacheRetries,
      avgTimeInCallCachingState = avgTimeInCallCaching,
      avgTimeInJobPreparation = avgTimeInJobPreparationState,
      avgTimeForFetchingAndCopyingCacheHit = avgTimeInCacheCopying
    )
  }
  else WorkflowMetrics(
    workflowId = workflowMetadata.id,
    workflowName = workflowMetadata.workflowName,
    workflowStartedAfter = workflowStartedAfterTime,
    workflowRunningTime = workflowRunningTimeDuration,
    totalJobsPerRootWf = totalJobsByRootWf,
    avgCacheCopyRetries = -1,
    avgTimeInCallCachingState = Duration.ZERO,
    avgTimeInJobPreparation = Duration.ZERO,
    avgTimeForFetchingAndCopyingCacheHit = Duration.ZERO
  )
}
