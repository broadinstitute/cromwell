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
  val bucket = params.reporterConfig.as[String]("report-bucket")
  lazy val cromwellVersion = CentaurCromwellClient.version.unsafeRunSync().cromwell

  /** A description of where the reporter is sending the errors. */
  override def destination = "GCS bucket"

  /** Send a report of a centaur failure. */
  override def logCentaurFailure(testEnvironment: TestEnvironment, ciEnvironment: CiEnvironment, centaurTestException: CentaurTestException)(implicit executionContext: ExecutionContext) = {
    centaurTestException.metadataJsonOption.map(pushJsonToGcs(_, centaurTestException.testName)).getOrElse(IO.unit)
  }

  override def logSuccessfulRun(testName: String, submitResponse: SubmitWorkflowResponse) = for {
    metadata <- CentaurCromwellClient.metadata(submitResponse.submittedWorkflow.id)
    _ <- pushJsonToGcs(metadata.value, testName)
  } yield ()

  private def pushJsonToGcs(json: String, testName: String) = IO {
    val path = makePath(testName)
    logger.info(s"Reporting metadata to gs://$bucket/$path")

    storage.create(
      BlobInfo.newBuilder(bucket, path)
        .setContentType("application/json")
        .build(),
      json.toArray.map(_.toByte)
    )
  }.void
  
  private def makePath(testName: String) = s"$testName/$cromwellVersion/metadata.json"
}
