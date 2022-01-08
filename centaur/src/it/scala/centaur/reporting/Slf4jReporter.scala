package centaur.reporting

import cats.effect.IO
import centaur.test.CentaurTestException
import com.typesafe.scalalogging.StrictLogging
import org.testcontainers.shaded.org.apache.commons.lang.exception.ExceptionUtils

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

  override def logFailure(testEnvironment: TestEnvironment,
                          ciEnvironment: CiEnvironment,
                          throwable: Throwable)
                         (implicit executionContext: ExecutionContext): IO[Unit] = {
    IO {

      val errorMessage = throwable match {
        case centaurTestException: CentaurTestException =>
          centaurTestException.workflowIdOption.map("with workflow id '" + _ + "' ").getOrElse("")
        case nonCentaurException =>
          s"with unexpected non-Centaur exception $nonCentaurException"
      }

      val message =
        s"Test '${testEnvironment.name}' " +
          s"failed on attempt ${testEnvironment.attempt + 1} " +
          s"of ${testEnvironment.retries + 1} " +
          errorMessage

      // Only log fully on the final attempt. Otherwise log a shortened version
      if (testEnvironment.attempt >= testEnvironment.retries) {
        logger.error(message, throwable)
      } else {
        val messageWithShortExceptionContext = message + " (" + ExceptionUtils.getMessage(throwable).replace("\n", " ").take(150) + "[...])"
        logger.warn(messageWithShortExceptionContext)
      }
    }
  }
}
