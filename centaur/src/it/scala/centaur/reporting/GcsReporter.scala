package centaur.reporting
import cats.effect.IO
import centaur.api.CentaurCromwellClient
import centaur.test.CentaurTestException
import centaur.test.submit.SubmitWorkflowResponse
import com.google.cloud.storage._
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext

class GcsReporter(override val params: ErrorReporterParams) extends ErrorReporter with SuccessReporter with StrictLogging {
  val storage = StorageOptions.getDefaultInstance.getService
  val reportBucket = params.reporterConfig.as[String]("report-bucket")
  val reportPath = params.reporterConfig.as[String]("report-path")

  /** A description of where the reporter is sending the errors. */
  override def destination = "GCS bucket"

  /**
    * In this ErrorReporter implementation this method will save information about exceptions of type
    * CentaurTestException to GCS. Exceptions of other types will be ignored.
    */
  override def logFailure(testEnvironment: TestEnvironment,
                          ciEnvironment: CiEnvironment,
                          throwable: Throwable)
                         (implicit executionContext: ExecutionContext): IO[Unit] = {
    throwable match {
      case centaurTestException: CentaurTestException =>
        logger.info(s"Reporting failed metadata to gs://$reportBucket/$reportPath")
        centaurTestException.metadataJsonOption.map(pushJsonToGcs).getOrElse(IO.unit)
      case _ =>
        IO.unit // this ErrorReporter only supports exceptions of CentaurTestException type
    }
  }

  override def logSuccessfulRun(submitResponse: SubmitWorkflowResponse): IO[Unit] = {
    logger.info(s"Reporting successful metadata to gs://$reportBucket/$reportPath")
    for {
      metadata <- CentaurCromwellClient.metadataWithId(submitResponse.submittedWorkflow.id)
      _ <- pushJsonToGcs(metadata.value)
    } yield ()
  }

  private def pushJsonToGcs(json: String) = IO {
    storage.create(
      BlobInfo.newBuilder(reportBucket, reportPath)
        .setContentType("application/json")
        .build(),
      json.toArray.map(_.toByte)
    )
  }.void
}
