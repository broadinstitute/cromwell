package cromwell.backend.impl.jes

import java.time.OffsetDateTime

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
    def eventList: Seq[EventStartTime]
    def machineType: Option[String]
    def zone: Option[String]
    def instanceName: Option[String]
  }

  case class Success(eventList: Seq[EventStartTime], machineType: Option[String], zone: Option[String], instanceName: Option[String]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  final case class Failed(errorCode: Int, errorMessage: Option[String], eventList: Seq[EventStartTime], machineType: Option[String], zone: Option[String], instanceName: Option[String])
    extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }
}

// An event with a startTime timestamp
case class EventStartTime(name: String, offsetDateTime: OffsetDateTime)
