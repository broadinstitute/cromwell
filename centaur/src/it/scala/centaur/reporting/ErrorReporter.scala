package centaur.reporting

import cats.effect.IO
import centaur.test.CentaurTestException

/**
  * Reports errors during testing.
  */
trait ErrorReporter {
  /** The name of this reporter. */
  def name: String

  /** A description of where the reporter is sending the errors. */
  def destination: String

  /** Send a report of a centaur failure. */
  def logCentaurFailure(testEnvironment: TestEnvironment,
                        ciEnvironment: CiEnvironment,
                        centaurTestException: CentaurTestException): IO[Unit]
}
