package cromwell.backend.io

import cromwell.core.{CallContext, StandardPaths}
import cromwell.core.path.Path

object JobPathsSpecHelper {
  implicit class EnhancedJobPaths(val jobPaths: JobPaths) extends AnyVal {
    def stdout: Path = jobPaths.standardPaths.output
    def stderr: Path = jobPaths.standardPaths.error
  }

  implicit class EnhancedCallContext(val callContext: CallContext) extends AnyVal {
    def stdout: String = callContext.standardPaths.output.pathAsString
    def stderr: String = callContext.standardPaths.error.pathAsString
  }

  val DummyStandardPaths: StandardPaths = null
}
