package cromwell.engine.backend

import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.GoogleAuthMode
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions

import scala.util.Try

object EnhancedWorkflowOptions {

  implicit class GoogleAuthWorkflowOptions(val workflowOptions: WorkflowOptions) extends AnyVal {
    def toGoogleAuthOptions: GoogleAuthMode.GoogleAuthOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }
  }
}
