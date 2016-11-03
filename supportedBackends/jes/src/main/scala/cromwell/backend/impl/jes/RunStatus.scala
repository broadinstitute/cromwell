package cromwell.backend.impl.jes

import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.{GlobComplete, GlobVerification}
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
  case class AwaitingGlobConsistency(globChecks: Seq[GlobVerification], eventList: Seq[ExecutionEvent], machineType: Option[String], zone: Option[String], instanceName: Option[String]) extends RunStatus {
    override def toString = s"AwaitingGlobConsistency(${globChecks.mkString(", ")})"
  }

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
    def machineType: Option[String]
    def zone: Option[String]
    def instanceName: Option[String]
  }

  case class Success(precalculatedGlobs: Seq[GlobComplete], eventList: Seq[ExecutionEvent], machineType: Option[String], zone: Option[String], instanceName: Option[String]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  final case class Failed(errorCode: Int, errorMessage: List[String], eventList: Seq[ExecutionEvent], machineType: Option[String], zone: Option[String], instanceName: Option[String])
    extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }
}
