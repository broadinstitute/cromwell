package centaur.reporting
import cats.effect.IO
import cats.syntax.functor._
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

  /** Send a report of a centaur failure. */
  override def logCentaurFailure(testEnvironment: TestEnvironment, ciEnvironment: CiEnvironment, centaurTestException: CentaurTestException)(implicit executionContext: ExecutionContext) = {
    logger.info(s"Reporting failed metadata to gs://$reportBucket/$reportPath")
    centaurTestException.metadataJsonOption.map(pushJsonToGcs).getOrElse(IO.unit)
  }

  override def logSuccessfulRun(submitResponse: SubmitWorkflowResponse) = {
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
