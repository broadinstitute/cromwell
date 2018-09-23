#!/usr/bin/env amm

@
import $ivy.`io.spray::spray-json:1.3.4`
import $ivy.`com.github.pathikrit::better-files:2.17.1`
import $ivy.`joda-time:joda-time:2.10`
import $ivy.`org.apache.commons:commons-lang3:3.8`
import $ivy.`com.github.tototoshi::scala-csv:1.3.5`

import spray.json._
import spray.json.DefaultJsonProtocol
import better.files.File
import java.time.{Duration, OffsetDateTime}
import java.time.format.DateTimeFormatter

import org.apache.commons.lang3.time.DurationFormatUtils
import com.github.tototoshi.csv._

import scala.collection.immutable.ListMap

//case class CallCachingHitFailuresMsg(causedBy: Option[Seq[CallCachingHitFailuresMsg]],
//                                     message: String)

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
                           totalJobsPerRootWf: Int)

case class ScatterWidthMetrics(taskName: String,
                               width: Int)

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
def extractMetricsFromMetadata(): Unit = {
  import MetadataJsonProtocol._

  def diffBetweenDateTimeInHumanReadableFormat(date1: OffsetDateTime, date2: OffsetDateTime): String = {
    val durationBetweenDates = Duration.between(date1, date2)
    DurationFormatUtils.formatDurationWords(durationBetweenDates.toMillis, true, true)
  }

  def convertCaseClassToMap(cc: AnyRef) =
    (ListMap[String, Any]() /: cc.getClass.getDeclaredFields) {
      (a, f) =>
        f.setAccessible(true)
        a + (f.getName -> f.get(cc))
    }


  val dateTimeFormatterPattern = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:SS")

  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/small_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/hello_world_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/subworkflow_hello_world_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/subworkflow_with_scatter_metadata.json")
  //  val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/cc_after_metadata.json")
  //    val metadataFile1 = File("/Users/sshah/Documents/perf_metadata_compare/large_scatter_with_multiple_calls_metadata.json")

  val metadataFile1Content = metadataFile1.contentAsString
  val workflowMetadata = metadataFile1Content.parseJson.convertTo[Metadata]

  println(s"************ WORKFLOW RUNNING TIMES **************")
  println(s"Submission time: ${workflowMetadata.submission}")
  println(s"Start time: ${workflowMetadata.start}")
  println(s"End time: ${workflowMetadata.end}")

  val workflowStartedAfterTime = diffBetweenDateTimeInHumanReadableFormat(workflowMetadata.submission, workflowMetadata.start)
  val workflowRunningTime = diffBetweenDateTimeInHumanReadableFormat(workflowMetadata.start, workflowMetadata.end)

  println(s"\nWorkflow was started after $workflowStartedAfterTime")
  println(s"Workflow ran for $workflowRunningTime")

  //DONE START- REMOVE IT
  println(s"*********** CACHE COPY TRIES BEFORE FINDING A SUCCESSFUL HIT (PER JOB) *************")
  workflowMetadata.calls.foreach(calls => {
    calls.foreach(jobsForCall => {
      val jobName = jobsForCall._1
      println(s"Job Name: $jobName")
      jobsForCall._2.foreach(job => {
        job.callCaching.map(callCacheObject => {
          val hitFailureMap = callCacheObject.hitFailures.getOrElse(Map.empty)
          println(s"\t-Shard: ${job.shardIndex} retires: ${hitFailureMap.size}")
        })

      })
    })
  })
  //DONE END- REMOVE IT

  //DONE START- BUT KEEP THIS
  var totalJobsByRootWf: Int = 0
  workflowMetadata.calls.foreach(callMap => {
    callMap.foreach(task => totalJobsByRootWf += task._2.size)
  })
  //DONE END- BUT KEEP THIS

  //DONE START- BUT KEEP THIS
  println(s"********** SCATTER WIDTH AT EXPANSION TIME **********")
  val scatterWidthMetricsOption = workflowMetadata.calls.map(callMap => callMap.map(task => {
    if(task._2.head.shardIndex != -1){
      println(s"Job ${task._1} scattered ${task._2.size} wide")
      ScatterWidthMetrics(task._1, task._2.size)
    }
    else {
      println(s"Job ${task._1} didn't have a scatter")
      ScatterWidthMetrics(task._1, -1)
    }
  }))
  //DONE END- BUT KEEP THIS

  //DONE START- REMOVE IT
  println(s"************* TIME BETWEEN SUBMISSION & CC HIT BEING SUCCESSFUL *********")
  val iWantTheseStates = List("CheckingCallCache", "FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")

  workflowMetadata.calls.foreach(c => c.foreach(s=> {
    println(s"Job Name: ${s._1}")
    s._2.foreach(call => {
      val eventsRelatedToCC = call.executionEvents.filter(e => iWantTheseStates.exists(s => s.equalsIgnoreCase(e.description)))

      if(eventsRelatedToCC.nonEmpty){
        val sortedEventsByStartTime = eventsRelatedToCC.sortBy(a => a.startTime)

        val timeInCallCachingState = diffBetweenDateTimeInHumanReadableFormat(sortedEventsByStartTime(0).startTime, sortedEventsByStartTime(sortedEventsByStartTime.size-1).endTime)

        println(s"\t-Shard ${call.shardIndex} spent $timeInCallCachingState in Call Caching state")
      }else
        println(s"\t-Shard ${call.shardIndex} did not do Call Caching")
    })}))
  //DONE END- REMOVE IT


  //NEW CODE- KEEP IT
  val callCachingEventStates = List("CheckingCallCache", "FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")
  val jobPreparationEventStates = List("Pending", "RequestingExecutionToken", "WaitingForValueStore", "PreparingJob", "CheckingJobStore")
  val cacheCopyingEventStates = List("FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")

  val taskLevelMetricsOption: Option[Iterable[TaskMetrics]] = workflowMetadata.calls.map(taskMap => taskMap.map(task => {
    val allJobMetrics: Seq[JobMetrics] = task._2.map(call => {
      //cache copy tries before finding a successful hit
      val hitFailuresMap = call.callCaching.map(callCachingObject => callCachingObject.hitFailures.getOrElse(Map.empty))

      val cacheRetries = if(hitFailuresMap.isDefined) hitFailuresMap.size else -1

      //time job spent in call caching states
      val eventsRelatedToCC = call.executionEvents.filter(event => callCachingEventStates.exists(state => state.equalsIgnoreCase(event.description)))

      val timeInCallCachingState = if (eventsRelatedToCC.nonEmpty) {
        val eventsSortedByStartTime = eventsRelatedToCC.sortBy(x => x.startTime)
        diffBetweenDateTimeInHumanReadableFormat(eventsSortedByStartTime.head.startTime, eventsSortedByStartTime.last.endTime)
      }
      else ""

      //time job spent in job preparation state i.e time between job submission and running it
      //(this doesn't consider the time it spent in call caching states
      val durationOfEventsInPreparationState = call.executionEvents
        .filter(event => jobPreparationEventStates.exists(state => state.equalsIgnoreCase(event.description)))
        .map(event => Duration.between(event.startTime, event.endTime))

      val timeInPreparationState = if(durationOfEventsInPreparationState.nonEmpty) {
        val totalTimeInPreparation = durationOfEventsInPreparationState.reduce(_ plus _)
        DurationFormatUtils.formatDurationWords(totalTimeInPreparation.toMillis, true, true)
      }
      else ""

      //time job spent in fetching and copying cache hit(s)
      val durationOfEventsInCacheCopyingState = call.executionEvents
        .filter(event => cacheCopyingEventStates.exists(state => state.equalsIgnoreCase(event.description)))
        .map(event => Duration.between(event.startTime, event.endTime))

      val timeInCopyingCache = if(durationOfEventsInCacheCopyingState.nonEmpty){
        val totalTimeCopyingCacheHits = durationOfEventsInCacheCopyingState.reduce(_ plus _)
        DurationFormatUtils.formatDurationWords(totalTimeCopyingCacheHits.toMillis, true, true)
      }
      else ""


      JobMetrics(
        shard = call.shardIndex,
        cacheCopyRetries = cacheRetries,
        timeInCallCachingState = timeInCallCachingState,
        timeFromSubmissionToRunning = timeInPreparationState,
        timeForFetchingCopyingCacheHit = timeInCopyingCache
      )
    })

    TaskMetrics(
      taskName = task._1,
      jobMetrics = allJobMetrics
    )
  }))
  //NEW CODE- KEEP IT



  //DONE START- REMOVE IT
  println(s"*********** TIME BETWEEN SUBMISSION OF JOB AND RUNNING IT ***********")
  val preparingJobEventsList = List("Pending", "RequestingExecutionToken", "WaitingForValueStore", "PreparingJob", "CheckingJobStore")
  workflowMetadata.calls.foreach(callArray => callArray.foreach(s => {
    println(s"Job Name: ${s._1}")
    s._2.foreach(call => {
      val eventsForThisJob: Seq[Duration] = call.executionEvents
        .filter(e => preparingJobEventsList.exists(s => s.equalsIgnoreCase(e.description)))
        .map(e => Duration.between(e.startTime, e.endTime))

      if(eventsForThisJob.nonEmpty){
        val totalTimeSpentBeforeRunning = eventsForThisJob.reduce(_ plus _)
        val textRep = DurationFormatUtils.formatDurationWords(totalTimeSpentBeforeRunning.toMillis, true, true)

        println(s"\t-Shard ${call.shardIndex}: $textRep")
      }
      else {
        println(s"\t-Shard ${call.shardIndex} somehow didn't spend time preparing for the job")
      }
    })}))
  //DONE END- REMOVE IT

  //DONE START- REMOVE IT
  println(s"*********** TIME FETCHING AND COPYING CAHE HIT(s) FOR A JOB **********")
  val copyingEvents = List("FetchingCachedOutputsFromDatabase", "BackendIsCopyingCachedOutputs")
  workflowMetadata.calls.foreach(callArray => callArray.foreach(s => {
    println(s"Job Name: ${s._1}")
    s._2.foreach(call => {
      val eventsForThisJob: Seq[Duration] = call.executionEvents
        .filter(e => copyingEvents.exists(s => s.equalsIgnoreCase(e.description)))
        .map(e => Duration.between(e.startTime, e.endTime))

      if(eventsForThisJob.nonEmpty){
        val totalTimeCopyingCacheHits = eventsForThisJob.reduce(_ plus _)
        val textRep = DurationFormatUtils.formatDurationWords(totalTimeCopyingCacheHits.toMillis, true, true)

        println(s"\t-Shard ${call.shardIndex}: $textRep")
      }
      else {
        println(s"\t-Shard ${call.shardIndex} didn't spend time fetching & copying cache hits!")
      }
    })}))
  //DONE END- REMOVE IT


  //BELOW EVERYTHING IS TO BE KEPT
  val workflowMetricsObject = WorkflowMetrics(
    workflowId = workflowMetadata.id,
    workflowName = workflowMetadata.workflowName,
    workflowStartedAfter = workflowStartedAfterTime,
    workflowRunningTime = workflowRunningTime,
    totalJobsPerRootWf = totalJobsByRootWf
  )



  //------------------------------------------------------------------------------------------------------------------
  //WRITE TO CSV

  val outputFile = new java.io.File("/Users/sshah/Documents/perf_metadata_compare/small_metadata_metrics.csv")
  val csvWriter = CSVWriter.open(outputFile)

  val workflowMetricsFormat = convertCaseClassToMap(workflowMetricsObject).map(row => List(row._1, row._2.toString)).toSeq
  csvWriter.writeAll(workflowMetricsFormat)

  csvWriter.writeRow("**********************")

  scatterWidthMetricsOption.foreach( scatterWidthMetrics => {
    csvWriter.writeRow(List("Scatter Width at expansion time"))
    val scatterMapOfKeyValues = scatterWidthMetrics.map(convertCaseClassToMap(_))
    val values = scatterMapOfKeyValues.map(_.values.toSeq)
    val header = scatterMapOfKeyValues.head.keys.toSeq

    csvWriter.writeAll(Seq(header) ++ values)
  })

  csvWriter.writeRow("**********************")




}
