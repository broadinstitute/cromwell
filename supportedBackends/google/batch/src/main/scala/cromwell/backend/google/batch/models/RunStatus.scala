package cromwell.backend.google.batch.models

import com.google.cloud.batch.v1.JobStatus
import cromwell.core.ExecutionEvent
import org.slf4j.{Logger, LoggerFactory}

sealed trait RunStatus

object RunStatus {

  private val log: Logger = LoggerFactory.getLogger(RunStatus.toString)

  def fromJobStatus(status: JobStatus.State): RunStatus = status match {
    case JobStatus.State.QUEUED =>
      log.info("job queued")
      Running
    case JobStatus.State.SCHEDULED =>
      log.info("job scheduled")
      Running
    case JobStatus.State.RUNNING =>
      log.info("job running")
      Running
    case JobStatus.State.SUCCEEDED =>
      log.info("job scheduled")
      Succeeded(List(ExecutionEvent("complete in GCP Batch"))) // update to more specific
    case JobStatus.State.FAILED =>
      log.info("job failed")
      Failed(List.empty)
    case JobStatus.State.DELETION_IN_PROGRESS =>
      log.info("deletion in progress")
      DeletionInProgress
    case JobStatus.State.STATE_UNSPECIFIED =>
      log.info("state unspecified")
      StateUnspecified
    case JobStatus.State.UNRECOGNIZED =>
      log.info("state unrecognized")
      Unrecognized
    case _ =>
      log.info(s"job status not matched: $status")
      Running
  }

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus

  case object Running extends RunStatus
  case object DeletionInProgress extends RunStatus
  case object StateUnspecified extends RunStatus
  case object Unrecognized extends RunStatus

  case class Succeeded(override val eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
    override def toString = "Succeeded"
  }

  final case class Failed(override val eventList: Seq[ExecutionEvent]) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }
}
