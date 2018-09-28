package perf

import better.files.File
import io.circe
import io.circe.generic.auto._
import io.circe.parser._

// Do not remove this unused import
// If removed circe parser fails to find implicit decoder for OffsetDateTime, and parsing fails
import io.circe.java8.time.decodeOffsetDateTimeDefault

object CompareMetadata extends App {

  def parseMetadataFromFile(filePath: String): Either[circe.Error, Metadata] = {
//    val metadataFile = File("/Users/sshah/Documents/perf_metadata_compare/small_metadata.json")
    val metadataFile = File(filePath)
    val metadataFileContent = metadataFile.contentAsString

    decode[Metadata](metadataFileContent)
  }

  def computeMetricsFromMetadata(metadata: Metadata): Unit = {
    println(metadata.avgTimeInCallCachingState)
  }

  def compareMetadataMetrics(metadataOld: Metadata, metadataNew: Metadata): Unit = {
    if (metadataNew.workflowStartedAfter.compareTo(metadataOld.workflowStartedAfter.multipliedBy(1.1.toLong)) > 0) {
      Console.err.println(s"The time after which workflow was started of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.workflowStartedAfter}" +
        s"\nOld metric value: ${metadataOld.workflowStartedAfter}")
      System.exit(1)
    }

    if (metadataNew.workflowRunningTime.compareTo(metadataOld.workflowRunningTime.multipliedBy(1.1.toLong)) > 0) {
      Console.err.println(s"The workflow running time of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.workflowRunningTime}" +
        s"\nOld metric value: ${metadataOld.workflowRunningTime}")
      System.exit(1)
    }

    if (metadataNew.avgCacheRetries > 1.1 * metadataOld.avgCacheRetries) {
      Console.err.println(s"Avg. Cache Copy Retries of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.avgCacheRetries}" +
        s"\nOld metric value: ${metadataOld.avgCacheRetries}")
      System.exit(1)
    }

    if (metadataNew.avgTimeInCallCachingState.compareTo(metadataOld.avgTimeInCallCachingState.multipliedBy(1.1.toLong)) > 0) {
      Console.err.println(s"Avg. Call Caching time of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.avgTimeInCallCachingState}" +
        s"\nOld metric value: ${metadataOld.avgTimeInCallCachingState}")
      System.exit(1)
    }

    if (metadataNew.avgTimeInJobPreparation.compareTo(metadataOld.avgTimeInJobPreparation.multipliedBy(1.1.toLong)) > 0) {
      Console.err.println(s"Avg. time in Job Preparation of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.avgTimeInJobPreparation}" +
        s"\nOld metric value: ${metadataOld.avgTimeInJobPreparation}")
      System.exit(1)
    }

    if (metadataNew.avgTimeForFetchingAndCopyingCacheHit.compareTo(metadataOld.avgTimeForFetchingAndCopyingCacheHit.multipliedBy(1.1.toLong)) > 0) {
      Console.err.println(s"Avg. time in Fetching and Copying Cache Hits of new metadata is greater than 10% of old metadata." +
        s"\nNew metric value: ${metadataNew.avgTimeForFetchingAndCopyingCacheHit}" +
        s"\nOld metric value: ${metadataOld.avgTimeForFetchingAndCopyingCacheHit}")
      System.exit(1)
    }
  }

  def printParseErrorToConsoleAndExit(errorList: circe.Error*): Unit = {
    Console.err.println(s"Something went wrong while parsing metadata json. Error:")
    errorList.foreach(e => e.printStackTrace(Console.err))
    System.exit(1)
  }


  args.length match {
    case 2 => {
      val metadataOldEither = parseMetadataFromFile(args(0))
      val metadataNewEither = parseMetadataFromFile(args(1))
      (metadataOldEither, metadataNewEither) match {
        case (Right(metadataOld), Right(metadataNew)) => compareMetadataMetrics(metadataOld, metadataNew)
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
