package cromwell.backend.impl.jes

import java.time.OffsetDateTime

sealed trait RunStatus {
  import RunStatus._

  // Could be defined as false for Initializing and true otherwise, but this is more defensive.
  def isRunningOrComplete = this match {
    case _:Running | _: TerminalRunStatus => true
    case _ => false
  }

  def eventList: Seq[EventStartTime]
}

object RunStatus {
  case class Initializing(eventList: Seq[EventStartTime]) extends RunStatus
  case class Running(eventList: Seq[EventStartTime]) extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[EventStartTime]
  }

  case class Success(eventList: Seq[EventStartTime]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  final case class Failed(errorCode: Int, errorMessage: Option[String], eventList: Seq[EventStartTime])
    extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }
}

// An event with a startTime timestamp
case class EventStartTime(name: String, offsetDateTime: OffsetDateTime)
