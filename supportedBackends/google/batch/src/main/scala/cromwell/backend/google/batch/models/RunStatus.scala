package cromwell.backend.google.batch.models

import cromwell.core.ExecutionEvent
import io.grpc.Status
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

  sealed trait TerminalRunStatus extends RunStatus

  case class Success(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends TerminalRunStatus {
    override def toString = "Success"
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val errorCode: Status
  }

  final case class Failed(
    errorCode: Status,
    eventList: Seq[ExecutionEvent],
    instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Aborted(errorCode: Status, instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    override def eventList: Seq[ExecutionEvent] = List.empty
  }

  // TODO: Use this when detecting a preemption or remove it
  final case class Preempted(
    errorCode: Status,
    eventList: Seq[ExecutionEvent],
    instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends UnsuccessfulRunStatus {
    override def toString = "Preempted"
  }
}
