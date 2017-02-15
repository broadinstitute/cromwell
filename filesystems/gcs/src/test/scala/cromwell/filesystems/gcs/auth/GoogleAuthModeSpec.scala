package cromwell.filesystems.gcs.auth

import cromwell.core.WorkflowOptions
import org.scalatest.Assertions._

import scala.util.{Failure, Try}

object GoogleAuthModeSpec {
  def assumeHasApplicationDefaultCredentials(): Unit = {
    tryApplicationDefaultCredentials match {
      case Failure(exception) => cancel(exception.getMessage)
      case _ =>
    }
    ()
  }

  private lazy val tryApplicationDefaultCredentials: Try[Unit] = Try {
    val authMode = ApplicationDefaultMode("application-default")
    val workflowOptions = WorkflowOptions.empty
    authMode.credential(workflowOptions)
    ()
  }
}
