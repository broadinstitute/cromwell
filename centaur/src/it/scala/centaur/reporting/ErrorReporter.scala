package centaur.reporting

import cats.effect.IO
import centaur.test.CentaurTestException

import scala.concurrent.ExecutionContext

/**
  * Reports errors during testing.
  */
trait ErrorReporter {
  /** The various parameters for this reporter. */
  def params: ErrorReporterParams

  /** A description of where the reporter is sending the errors. */
  def destination: String

  /** Send a report of a centaur failure. */
  def logCentaurFailure(testEnvironment: TestEnvironment,
                        ciEnvironment: CiEnvironment,
                        centaurTestException: CentaurTestException)
                       (implicit executionContext: ExecutionContext): IO[Unit]
}
