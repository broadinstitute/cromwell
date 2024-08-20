package cromwell.backend.google.batch.models

import cromwell.core.ExecutionEvent

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
    val exitCode: Option[GcpBatchExitCode]
    val prettyPrintedError: String
  }

  final case class Failed(
    exitCode: Option[GcpBatchExitCode],
    eventList: Seq[ExecutionEvent]
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

  final case object Aborted extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    override val exitCode: Option[GcpBatchExitCode] = None

    override def eventList: Seq[ExecutionEvent] = List.empty

    override val prettyPrintedError: String = "The job was aborted"
  }
}
