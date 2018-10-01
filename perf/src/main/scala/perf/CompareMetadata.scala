package perf

import better.files.File
import io.circe
import io.circe.generic.auto._
import io.circe.parser._
import cats.implicits._
import java.time.Duration

import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr._

// Do not remove this unused import
// If removed, circe parser fails to find implicit decoder for OffsetDateTime, and parsing fails
import io.circe.java8.time.decodeOffsetDateTimeDefault

object CompareMetadata extends App {

  def parseMetadataFromFile(filePath: String): Either[circe.Error, Metadata] = {
    val metadataFile = File(filePath)
    val metadataFileContent = metadataFile.contentAsString

    decode[Metadata](metadataFileContent)
  }


  def displayComputedMetrics(metadata: Metadata, displayMsg: String): Unit = {
    println(displayMsg)

    println(s"Workflow started after: ${metadata.workflowStartedAfter}")
    println(s"Workflow Running time: ${metadata.workflowRunningTime}")
    println(s"Total jobs per root workflow: ${metadata.totalJobsPerRootWf}")
    println(s"Avg cache copy retries: ${metadata.avgCacheRetries}")
    println(s"Avg time job spent in Call Caching state: ${metadata.avgTimeInCallCachingState}")
    println(s"Avg time job spent in Job Preparation state: ${metadata.avgTimeInJobPreparation}")
    println(s"Avg time job spent in fetching and copying cache hit(s) state: ${metadata.avgTimeForFetchingAndCopyingCacheHit}")
  }

  def compareDurationMetrics(metadataOld: Metadata, metadataNew: Metadata, metricFunc: Metadata => Duration, metricName: String): ErrorOr[String] = {
    if(metricFunc(metadataNew).compareTo(metricFunc(metadataOld).multipliedBy(1.1.toLong)) > 0){
      (s"$metricName of new metadata is greater than 10% of old metadata. " +
        s"New metric value: ${metricFunc(metadataNew)}. " +
        s"Old metric value: ${metricFunc(metadataOld)}.").invalidNel[String]
    }
    else s"$metricName hasn't regressed yet".valid
  }


  def compareIntMetrics(metadataOld: Metadata, metadataNew: Metadata, metricFunc: Metadata => Int, metricName: String): ErrorOr[String] = {
    if(metricFunc(metadataNew) > 1.1 * metricFunc(metadataOld)){
      (s"$metricName of new metadata is greater than 10% of old metadata. " +
        s"New metric value: ${metricFunc(metadataNew)}. " +
        s"Old metric value: ${metricFunc(metadataOld)}.").invalidNel[String]
    }
    else s"$metricName hasn't regressed yet".valid
  }


  def compareMetadataMetrics(metadataOld: Metadata, metadataNew: Metadata): ErrorOr[List[String]] = {
    val durationFuncList: List[(Metadata => Duration, String)] = List((_.workflowStartedAfter, "workflowStartedAfter"),
      (_.workflowRunningTime, "workflowRunningTime"),
      (_.avgTimeInCallCachingState, "avgTimeInCallCachingState"),
      (_.avgTimeInJobPreparation, "avgTimeInJobPreparation"),
      (_.avgTimeForFetchingAndCopyingCacheHit, "avgTimeForFetchingAndCopyingCacheHit"))

    val intFuncList:  List[(Metadata => Int, String)] = List((_.totalJobsPerRootWf, "totalJobsPerRootWf"), (_.avgCacheRetries, "avgCacheRetries"))

    val durationMetricsComp = durationFuncList.map(x => compareDurationMetrics(metadataOld, metadataNew, x._1, x._2))
    val intMetricsComp = intFuncList.map(x => compareIntMetrics(metadataOld, metadataNew, x._1, x._2))
    (durationMetricsComp ::: intMetricsComp).sequence[ErrorOr, String]
  }


  def printParseErrorToConsoleAndExit(errorList: circe.Error*): Unit = {
    Console.err.println(s"Something went wrong while parsing metadata json. Error:")
    errorList.foreach(e => e.printStackTrace(Console.err))
    System.exit(1)
  }


  args.length match {
    case 2 => {
      //The args can be path to file on GCS, but currently it only considers local filesystem paths
      val metadataOldEither = parseMetadataFromFile(args(0))
      val metadataNewEither = parseMetadataFromFile(args(1))

      (metadataOldEither, metadataNewEither) match {
        case (Right(metadataOld), Right(metadataNew)) => {
          val metadataOldMsg = s"Metrics for metadata generated from ${args(0)}"
          val metadataNewMsg = s"Metrics for metadata generated from ${args(1)}"
          displayComputedMetrics(metadataOld, metadataOldMsg)
          displayComputedMetrics(metadataNew, metadataNewMsg)
          compareMetadataMetrics(metadataOld, metadataNew) match {
            case Valid(_) => Console.println("YAY!! Metrics from new metadata json haven't regressed!")
            case Invalid(listOfErrors) => {
              Console.err.println("Below metadata metrics have regressed:")
              Console.err.println(listOfErrors.toList.mkString("\n"))
              System.exit(1)
            }
          }
        }
        case (Right(_), Left(e)) => printParseErrorToConsoleAndExit(e)
        case (Left(e), Right(_)) => printParseErrorToConsoleAndExit(e)
        case (Left(e1), Left(e2)) => printParseErrorToConsoleAndExit(e1, e2)
      }
    }
    case _ => {
      Console.err.println("Please pass in 2 file paths!")
      System.exit(1)
    }
  }
}
