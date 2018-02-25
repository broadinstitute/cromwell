package cromwell.backend.impl.bcs

import cromwell.core.ExecutionEvent

import scala.util.{Failure, Success, Try}

sealed trait RunStatus {
  import RunStatus._
  val jobId: String
  val status: String

  def isTerminated: Boolean

  def isRunningOrComplete = this match {
    case _: Running | _: TerminalRunStatus => true
    case _ => false
  }

  override def toString = status
}

object RunStatus {
  final case class  Waiting(override val jobId: String) extends RunStatus {
    override val status = "Waiting"

    override def isTerminated: Boolean = false
  }

  final case class  Running(override val jobId: String) extends RunStatus {
    override val status = "Running"

    override def isTerminated: Boolean = false
  }

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]

    override def isTerminated: Boolean = true
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val errorMessage: Option[String]
    lazy val prettyPrintedError: String = errorMessage map { e => s" Message: $e" } getOrElse ""
  }

  final case class Finished(override val jobId: String, eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
    override val status = "Finished"
  }

  object UnsuccessfulRunStatus {
    def apply(jobId: String, status: String, errorMessage: Option[String], eventList: Seq[ExecutionEvent]): UnsuccessfulRunStatus = {
      if (status == "Stopped") {
        Stopped(jobId, errorMessage, eventList)
      } else {
        Failed(jobId, errorMessage, eventList)
      }
    }
  }

  final case class Failed(override val jobId: String,
                          errorMessage: Option[String],
                          eventList: Seq[ExecutionEvent]) extends UnsuccessfulRunStatus {
    override val status = "Failed"
  }

  final case class Stopped(override val jobId: String,
                           errorMessage: Option[String],
                           eventList: Seq[ExecutionEvent]) extends UnsuccessfulRunStatus {
      override val status = "Stopped"
  }
}

object RunStatusFactory {
    def getStatus(jobId: String, status: String, errorMessage: Option[String] = None, eventList: Option[Seq[ExecutionEvent]] = None): Try[RunStatus] = {
        import RunStatus._
        status match {
            case "Waiting" => Success(Waiting(jobId))
            case "Running" => Success(Running(jobId))
            case "Stopped" => Success(Stopped(jobId, errorMessage, eventList.getOrElse(Seq.empty)))
            case "Failed"  => Success(Failed(jobId, errorMessage, eventList.getOrElse(Seq.empty)))
            case "Finished" => Success(Finished(jobId, eventList.getOrElse(Seq.empty)))
            case _ => Failure(new RuntimeException(s"job {$jobId} turns to an invalid batchcompue status: {$status}"))
        }
    }
}
