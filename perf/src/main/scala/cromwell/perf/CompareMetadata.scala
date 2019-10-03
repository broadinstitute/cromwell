package cromwell.perf

import java.io.FileInputStream
import java.time.{Duration, OffsetDateTime}

import better.files.File
import cats.data.Validated.{Invalid, Valid}
import cats.implicits._
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.storage.StorageOptions
import com.typesafe.scalalogging.StrictLogging
import common.validation.ErrorOr._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser._

object CompareMetadata extends App with StrictLogging{
  private val REGRESSION_CONST = 1.1

  // If removed, IntelliJ (Oct '18) thinks the import isn't used.
  // Later the compiler fails to find a decoder for OffsetDateTime.
  lazy val anchorOffsetDateTime: Decoder[OffsetDateTime] = implicitly
  implicit lazy val decodeCall: Decoder[Call] = deriveDecoder
  implicit lazy val decodeCallCaching: Decoder[CallCaching] = deriveDecoder
  implicit lazy val decodeExecutionEvents: Decoder[ExecutionEvents] = deriveDecoder
  implicit lazy val decodeMetadata: Decoder[Metadata] = deriveDecoder

  def parseMetadataFromLocalFile(filePath: String): Either[Error, Metadata] = {
    val metadataFile = File(filePath)
    val metadataFileContent = metadataFile.contentAsString

    decode[Metadata](metadataFileContent)
  }


  def parseMetadataFromGcsFile(gcsUrl: String, pathToServiceAccount: String): Either[Error, Metadata] = {
    val gcsUrlArray = gcsUrl.replace("gs://", "").split("/", 2)
    val Array(gcsBucket, fileToBeLocalized) = gcsUrlArray

    val credentials = GoogleCredentials.fromStream(new FileInputStream(pathToServiceAccount))
    val storage = StorageOptions.newBuilder().setCredentials(credentials).build().getService
    val blob = storage.get(gcsBucket, fileToBeLocalized)
    val metadataFileContent = blob.getContent().map(_.toChar).mkString

    decode[Metadata](metadataFileContent)
  }


  def parseMetadata(inputFile: String, pathToServiceAccount: String): Either[Error, Metadata] = {
    if (inputFile.startsWith("gs://"))
      parseMetadataFromGcsFile(inputFile, pathToServiceAccount)
    else parseMetadataFromLocalFile(inputFile)
  }


  def displayComputedMetrics(metadata: Metadata, displayMsg: String): Unit = {
    logger.info(displayMsg)

    logger.info(s"Workflow started after: ${metadata.workflowStartedAfter}")
    logger.info(s"Workflow Running time: ${metadata.workflowRunningTime}")
    logger.info(s"Total jobs per root workflow: ${metadata.totalJobsPerRootWf}")
    logger.info(s"Avg cache copy retries: ${metadata.avgCacheRetries}")
    logger.info(s"Avg time job spent in Call Caching state: ${metadata.avgTimeInCallCachingState}")
    logger.info(s"Avg time job spent in just CheckingCallCache state: ${metadata.avgTimeInCheckingCallCacheState}")
    logger.info(s"Avg time job spent in Job Preparation state: ${metadata.avgTimeInJobPreparation}")
    logger.info(s"Avg time job spent in fetching and copying cache hit(s) state: ${metadata.avgTimeForFetchingAndCopyingCacheHit}")
  }


  /***
    * Compares the metrics in Metadata which have type 'Duration' (mostly the ones related to time metrics)
    */
  def compareDurationMetrics(metadataOld: Metadata, metadataNew: Metadata, metricFunc: Metadata => Duration, metricName: String): ErrorOr[String] = {
    if(metricFunc(metadataNew).toMillis > REGRESSION_CONST * metricFunc(metadataOld).toMillis){
      (s"$metricName of new metadata is greater than 10% of old metadata. " +
        s"New metric value: ${metricFunc(metadataNew)}. " +
        s"Old metric value: ${metricFunc(metadataOld)}.").invalidNel[String]
    }
    else s"$metricName hasn't regressed yet".valid
  }


  /***
    * Compares the metrics in Metadata which have type 'Int' (mostly the ones related to size metrics)
    */
  def compareIntMetrics(metadataOld: Metadata, metadataNew: Metadata, metricFunc: Metadata => Int, metricName: String): ErrorOr[String] = {
    if(metricFunc(metadataNew) > REGRESSION_CONST * metricFunc(metadataOld)){
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
      (_.avgTimeInCheckingCallCacheState, "avgTimeInCheckingCallCacheState"),
      (_.avgTimeInJobPreparation, "avgTimeInJobPreparation"),
      (_.avgTimeForFetchingAndCopyingCacheHit, "avgTimeForFetchingAndCopyingCacheHit"))

    val intFuncList:  List[(Metadata => Int, String)] = List((_.totalJobsPerRootWf, "totalJobsPerRootWf"), (_.avgCacheRetries, "avgCacheRetries"))

    val durationMetricsComp = durationFuncList.map(x => compareDurationMetrics(metadataOld, metadataNew, x._1, x._2))
    val intMetricsComp = intFuncList.map(x => compareIntMetrics(metadataOld, metadataNew, x._1, x._2))
    (durationMetricsComp ::: intMetricsComp).sequence[ErrorOr, String]
  }


  def printParseErrorToConsoleAndExit(metadataFile: String, error: Error, systemExit: Boolean): Unit = {
    logger.error(s"Something went wrong while parsing $metadataFile. Error: ${error.getLocalizedMessage}")
    if (systemExit) System.exit(1)
  }


  def generateAndCompareMetrics(metadataOldEither: Either[Error, Metadata], metadataNewEither: Either[Error, Metadata]): Unit = {
    (metadataOldEither, metadataNewEither) match {
      case (Right(metadataOld), Right(metadataNew)) =>
        val metadataOldMsg = s"Metrics for metadata generated from ${args(0)}"
        val metadataNewMsg = s"\nMetrics for metadata generated from ${args(1)}"
        displayComputedMetrics(metadataOld, metadataOldMsg)
        displayComputedMetrics(metadataNew, metadataNewMsg)
        compareMetadataMetrics(metadataOld, metadataNew) match {
          case Valid(_) => logger.info("\nYAY!! Metrics from new metadata json haven't regressed!")
          case Invalid(listOfErrors) =>
            logger.error("\nBelow metadata metrics have regressed:")
            logger.error(listOfErrors.toList.mkString("\n"))
            System.exit(1)
        }
      case (Right(_), Left(e)) => printParseErrorToConsoleAndExit(args(1), e, systemExit = true)
      case (Left(e), Right(_)) => printParseErrorToConsoleAndExit(args(0), e, systemExit = true)
      case (Left(e1), Left(e2)) =>
        printParseErrorToConsoleAndExit(args(0), e1, systemExit = false)
        printParseErrorToConsoleAndExit(args(1), e2, systemExit = true)
    }
  }


  args.length match {
    case 2 =>
      if (args(0).startsWith("gs://") || args(1).startsWith("gs://")) {
        logger.error("Path to service account is needed to download GCS file. Please pass it as 3rd argument.")
        System.exit(1)
      }
      else generateAndCompareMetrics(parseMetadataFromLocalFile(args(0)), parseMetadataFromLocalFile(args(1)))
    case 3 => generateAndCompareMetrics(parseMetadata(args(0), args(2)), parseMetadata(args(1), args(2)))
    case _ =>
      logger.error("Please pass in 2 file paths!")
      System.exit(1)
  }
}
