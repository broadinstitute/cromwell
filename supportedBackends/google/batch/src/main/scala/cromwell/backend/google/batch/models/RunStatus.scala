package cromwell.backend.google.batch.models

import cromwell.core.ExecutionEvent
import io.grpc.Status

sealed trait RunStatus

object RunStatus {

  case object Initializing extends RunStatus

  case object AwaitingCloudQuota extends RunStatus

  case object Running extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
  }

  case class Success(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val errorCode: Status
  }

  final case class Failed(
    errorCode: Status,
    eventList: Seq[ExecutionEvent]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Aborted(errorCode: Status) extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    override def eventList: Seq[ExecutionEvent] = List.empty
  }

  final case class Preempted(
    errorCode: Status,
    eventList: Seq[ExecutionEvent]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Preempted"
  }
}
