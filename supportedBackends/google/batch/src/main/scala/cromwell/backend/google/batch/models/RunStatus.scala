package cromwell.backend.google.batch.models

import cromwell.core.ExecutionEvent
import cromwell.services.cost.InstantiatedVmInfo

sealed trait RunStatus {
  def eventList: Seq[ExecutionEvent]
  def toString: String

  val instantiatedVmInfo: Option[InstantiatedVmInfo]
}

object RunStatus {

  case class Initializing(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends RunStatus { override def toString = "Initializing" }
  case class AwaitingCloudQuota(eventList: Seq[ExecutionEvent],
                                instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends RunStatus {
    override def toString = "AwaitingCloudQuota"
  }

  case class Running(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends RunStatus { override def toString = "Running" }

  case class Aborting(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty) extends RunStatus { override def toString = "Aborting" }

  sealed trait TerminalRunStatus extends RunStatus

  case class Success(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends TerminalRunStatus {
    override def toString = "Success"
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus

  final case class Failed(
    errorCode: GcpBatchExitCode,
    eventList: Seq[ExecutionEvent],
    instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Aborted(instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    override def eventList: Seq[ExecutionEvent] = List.empty
  }
}
