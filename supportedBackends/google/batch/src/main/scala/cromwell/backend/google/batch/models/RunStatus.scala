package cromwell.backend.google.batch.models

import cromwell.core.ExecutionEvent
import io.grpc.Status

import scala.util.Try

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
    val errorMessages: List[String]
    lazy val prettyPrintedError: String = errorMessages.mkString("\n")
    val errorCode: Status

    /**
     * If one exists, the JES error code (not the google RPC) (extracted from the error message)
     */
    val jesCode: Option[Int]
  }

  object UnsuccessfulRunStatus {
    def apply(errorCode: Status,
              errorMessage: Option[String],
              eventList: Seq[ExecutionEvent]
    ): UnsuccessfulRunStatus = {
      val jesCode: Option[Int] = errorMessage.flatMap { em =>
        Try(em.substring(0, em.indexOf(':')).toInt).toOption
      }
      Failed(errorCode, jesCode, errorMessage.toList, eventList)
    }
  }

  final case class Failed(
    errorCode: Status,
    jesCode: Option[Int],
    errorMessages: List[String],
    eventList: Seq[ExecutionEvent]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Aborted(errorCode: Status) extends UnsuccessfulRunStatus {
    override def toString = "Aborted"

    // In Batch, a workflow can't be aborted but its deleted instead, due to this, we don't have any info about the job
    override val errorMessages: List[String] = List.empty
    override val jesCode: Option[Int] = None
    override def eventList: Seq[ExecutionEvent] = List.empty
  }

  final case class Preempted(
    errorCode: Status,
    jesCode: Option[Int],
    errorMessages: List[String],
    eventList: Seq[ExecutionEvent]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Preempted"
  }
}
