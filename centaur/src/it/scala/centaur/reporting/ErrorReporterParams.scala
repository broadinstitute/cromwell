package centaur.reporting

import com.typesafe.config.Config

/**
  * Collects all of the parameters to pass to a new ErrorReporter.
  */
case class ErrorReporterParams(
  name: String,
  rootConfig: Config,
  reporterConfig: Config,
  database: ErrorReporterCromwellDatabase
)
