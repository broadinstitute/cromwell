package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.JobStatus
import cromwell.core.ExecutionEvent
import org.slf4j.{Logger, LoggerFactory}

sealed trait RunStatus {
    //import RunStatus._

    //def isTerminal: Boolean

    /*
    def isRunningOrComplete = this match {
        case Running | _: TerminalRunStatus => true
        case _ => false
    }
    */
}

object RunStatus {

    private val log: Logger = LoggerFactory.getLogger(RunStatus.toString)


    //def fromJobStatus(status: JobStatus.State,  eventList: Seq[ExecutionEvent] = Seq.empty): Try[RunStatus] = {
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
            Succeeded(List(ExecutionEvent("complete in GCP Batch"))) //update to more specific
        case JobStatus.State.FAILED =>
            log.info("job failed")
            Failed
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
            log.info("job status not matched")
            Running
    }

    sealed trait TerminalRunStatus extends RunStatus {
        def eventList: Seq[ExecutionEvent]
    }

    case object Initializing extends RunStatus {
        //def isTerminal=false
    }
    case object Running extends RunStatus {
        //def isTerminal=false
    }
    case object Succeeded extends TerminalRunStatus {
        override def eventList: Seq[ExecutionEvent] = List.empty
        //def isTerminal=true
    }

    case object Failed extends TerminalRunStatus {

        override def eventList: Seq[ExecutionEvent] = List.empty
    }

    case object DeletionInProgress extends RunStatus

    case object StateUnspecified extends RunStatus
    case object Unrecognized extends RunStatus

    //case class Succeeded(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {

    case class Succeeded(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
        //machineType: Option[String],
        //zone: Option[String])
        //instanceName: Option[String]) extends TerminalRunStatus {
        override def toString = "Succeeded"
        //def isTerminal = true
    }

    final case class Failed(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
        override def toString = "Failed"
    }

    final case class Cancelled() {
        override def toString = "Cancelled"
    }
    final case class Preempted() {
        override def toString = "Preempted"
    }

    final case class DeletionInProgress() {
        override def toString = "Deletion in Progress"
    }

    final case class StateUnspecified() {
        override def toString = "State Unspecified"
    }

    final case class Unrecognized() {
        override def toString = "Unrecognized"
    }

}
