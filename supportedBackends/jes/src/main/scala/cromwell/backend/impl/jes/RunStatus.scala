package cromwell.backend.impl.jes

import cromwell.core.ExecutionEvent

import scala.util.Try

sealed trait RunStatus {
  import RunStatus._

  // Could be defined as false for Initializing and true otherwise, but this is more defensive.
  def isRunningOrComplete = this match {
    case Running | _: TerminalRunStatus => true
    case _ => false
  }
}

object RunStatus {
  case object Initializing extends RunStatus
  case object Running extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
    def machineType: Option[String]
    def zone: Option[String]
    def instanceName: Option[String]
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val errorMessage: Option[String]
    lazy val prettyPrintedError: String = errorMessage map { e => s" Message: $e" } getOrElse ""
    val errorCode: Int

    /**
      * If one exists, the JES error code (not the google RPC) (extracted from the error message)
      */
    val jesCode: Option[Int]
  }

  case class Success(eventList: Seq[ExecutionEvent],
                     machineType: Option[String],
                     zone: Option[String],
                     instanceName: Option[String]) extends TerminalRunStatus {
    override def toString = "Success"
  }

  object UnsuccessfulRunStatus {
    def apply(errorCode: Int,
              errorMessage: Option[String],
              eventList: Seq[ExecutionEvent],
              machineType: Option[String],
              zone: Option[String],
              instanceName: Option[String]): UnsuccessfulRunStatus = {
      val jesCode: Option[Int] = errorMessage flatMap { em => Try(em.substring(0, em.indexOf(':')).toInt).toOption }
      if (errorCode == JesAsyncBackendJobExecutionActor.GoogleAbortedRpc && jesCode.contains(JesAsyncBackendJobExecutionActor.JesPreemption)) {
        Preempted(errorCode, jesCode, errorMessage, eventList, machineType, zone, instanceName)
      } else {
        Failed(errorCode, jesCode, errorMessage, eventList, machineType, zone, instanceName)
      }
    }
  }

  final case class Failed(errorCode: Int,
                          jesCode: Option[Int],
                          errorMessage: Option[String],
                          eventList: Seq[ExecutionEvent],
                          machineType: Option[String],
                          zone: Option[String],
                          instanceName: Option[String]) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Preempted(errorCode: Int,
                             jesCode: Option[Int],
                             errorMessage: Option[String],
                             eventList: Seq[ExecutionEvent],
                             machineType: Option[String],
                             zone: Option[String],
                             instanceName: Option[String]) extends UnsuccessfulRunStatus {
    override def toString = "Preempted"
  }
}
