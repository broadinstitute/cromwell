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

  sealed trait TerminalRunStatus extends RunStatus

  case class Success(eventList: Seq[ExecutionEvent], instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty)
      extends TerminalRunStatus {
    override def toString = "Success"
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val exitCode: Option[GcpBatchExitCode]
    val prettyPrintedError: String
  }

  final case class Failed(
    exitCode: Option[GcpBatchExitCode],
    eventList: Seq[ExecutionEvent],
    instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"

    override val prettyPrintedError: String =
      exitCode match {
        case Some(code) =>
          code match {
            case GcpBatchExitCode.VMPreemption => "A Spot VM for the job was preempted during run time"
            case GcpBatchExitCode.VMReportingTimeout =>
              "There was a timeout in the backend that caused a VM for the job to no longer receive updates"
            case GcpBatchExitCode.VMRebootedDuringExecution => "A VM for the job unexpectedly rebooted during run time"
            case GcpBatchExitCode.VMAndTaskAreUnresponsive =>
              "A task reached the unresponsive time limit and cannot be cancelled"
            case GcpBatchExitCode.TaskRunsOverMaximumRuntime =>
              "A task's run time exceeded the time limit specified in the maxRunDuration, or, a runnable's run time exceeded the time limit specified in the timeout"
            case GcpBatchExitCode.VMRecreatedDuringExecution =>
              "A VM for a job is unexpectedly recreated during run time"
          }
        case None =>
          // Take the last event as that is more likely to be indicative of what killed the job than the first event.
          eventList.lastOption
            .map(_.name)
            .getOrElse(
              "The job has failed but the exit code couldn't be derived, there isn't an event message either, please review the logs and report a bug"
            )
      }
  }

  final case class Aborted(eventList: Seq[ExecutionEvent],
                           instantiatedVmInfo: Option[InstantiatedVmInfo] = Option.empty
  ) extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    override val exitCode: Option[GcpBatchExitCode] = None

    override val prettyPrintedError: String = "The job was aborted"
  }
}
