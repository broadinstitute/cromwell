package centaur.reporting

import cats.effect.IO
import centaur.test.CentaurTestException
import com.typesafe.scalalogging.StrictLogging

import scala.concurrent.ExecutionContext

/**
  * An error reporter that only prints to slf4j. These errors have a very high chance of being missed and ignored.
  *
  * Useful as a backup in cases where another reporter is not available, for example in external PRs where secure
  * environment variables are not available.
  */
class Slf4jReporter(override val params: ErrorReporterParams)
  extends ErrorReporter with StrictLogging {

  override lazy val destination: String = "error"

  override def logCentaurFailure(testEnvironment: TestEnvironment,
                                 ciEnvironment: CiEnvironment,
                                 centaurTestException: CentaurTestException)
                                (implicit executionContext: ExecutionContext): IO[Unit] = {
    IO {
      val message =
        s"Test '${testEnvironment.name}' " +
          centaurTestException.workflowIdOption.map("with workflow id '" + _ + "' ").getOrElse("") +
          s"failed on attempt ${testEnvironment.attempt + 1} " +
          s"of ${testEnvironment.retries + 1}"
      val metadata = centaurTestException.metadataJsonOption.map(" with metadata:\n" + _).getOrElse("")
      logger.error(message + metadata, centaurTestException)
    }
  }
}
