package centaur.reporting

import cats.effect.IO

import scala.concurrent.ExecutionContext

/**
  * Reports errors during testing.
  */
trait ErrorReporter {
  /** The various parameters for this reporter. */
  def params: ErrorReporterParams

  /** A description of where the reporter is sending the errors. */
  def destination: String

  /** Send a report of a failure. */
  def logFailure(testEnvironment: TestEnvironment,
                 ciEnvironment: CiEnvironment,
                 throwable: Throwable)
                (implicit executionContext: ExecutionContext): IO[Unit]
}
