package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.impl.jes.errors.JesError
import cromwell.core.ExecutionEvent

sealed trait RunStatus {
  import RunStatus._

  // Could be defined as false for Initializing and true otherwise, but this is more defensive.
  def isRunningOrComplete = this match {
    case Running | _: TerminalRunStatus => true
    case _ => false
  }
}

object RunStatus {
  case object Initializing extends RunStatus
  case object Running extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
    def machineType: Option[String]
    def zone: Option[String]
    def instanceName: Option[String]
  }

  case class Success(eventList: Seq[ExecutionEvent], machineType: Option[String], zone: Option[String], instanceName: Option[String]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  final case class Failed(errorCode: Int, errorMessage: Option[String], eventList: Seq[ExecutionEvent], machineType: Option[String], zone: Option[String], instanceName: Option[String])
    extends TerminalRunStatus {
    
    private def unknownError(jobTag: String) = {
      new RuntimeException(s"Task $jobTag failed: error code $errorCode. Message: $errorMessage")
    }
    
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
    
    def toFailure(jobTag: String, stderrPath: Option[Path]): Exception = {
      JesError.fromFailedStatus(this, jobTag, stderrPath) getOrElse unknownError(jobTag)
    }
  }
}
