package cromwell.backend.impl.local

import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.GoogleAuthMode
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions

import scala.util.Try

object LocalImplicits {

  // TODO how to harmonize this ? Already exists in JES but needs gcs dependency - so can't be in core or backend at the moment
  implicit class GoogleAuthWorkflowOptions(val workflowOptions: WorkflowOptions) extends AnyVal {
    def toGoogleAuthOptions: GoogleAuthMode.GoogleAuthOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }
  }

}
