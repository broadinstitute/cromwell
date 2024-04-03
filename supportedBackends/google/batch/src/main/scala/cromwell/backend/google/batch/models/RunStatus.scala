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
      val jesCode: Option[Int] = errorMessage flatMap { em => Try(em.substring(0, em.indexOf(':')).toInt).toOption }

      val unsuccessfulStatusBuilder = errorCode match {
        case Status.CANCELLED => Cancelled.apply _
        case _ => Failed.apply _
      }

      unsuccessfulStatusBuilder.apply(errorCode, jesCode, errorMessage.toList, eventList)
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

  /**
   * What Cromwell calls Aborted, PAPI calls Cancelled. This means the job was "cancelled" by the user
   */
  final case class Cancelled(
    errorCode: Status,
    jesCode: Option[Int],
    errorMessages: List[String],
    eventList: Seq[ExecutionEvent]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Cancelled"
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
