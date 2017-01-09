package cromwell.filesystems.gcs.auth

import cromwell.core.WorkflowOptions
import org.scalatest.Assertions._

object GoogleAuthModeSpec {
  def assumeHasApplicationDefaultCredentials(): Unit = {
    try {
      val authMode = ApplicationDefaultMode("application-default")
      val workflowOptions = WorkflowOptions.empty
      authMode.authCredentials(workflowOptions)
      authMode.credential(workflowOptions)
      ()
    } catch {
      case exception: Exception => cancel(exception.getMessage)
    }
  }
}
