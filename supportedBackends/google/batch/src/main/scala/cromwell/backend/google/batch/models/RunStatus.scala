package cromwell.backend.google.batch.models

import cromwell.backend.google.batch.actors.GcpBatchAsyncBackendJobExecutionActor
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

    def machineType: Option[String]

    def zone: Option[String]

    def instanceName: Option[String]
  }

  case class Success(eventList: Seq[ExecutionEvent],
                     machineType: Option[String],
                     zone: Option[String],
                     instanceName: Option[String]
  ) extends TerminalRunStatus {
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
              eventList: Seq[ExecutionEvent],
              machineType: Option[String],
              zone: Option[String],
              instanceName: Option[String],
              wasPreemptible: Boolean
    ): UnsuccessfulRunStatus = {
      val jesCode: Option[Int] = errorMessage flatMap { em => Try(em.substring(0, em.indexOf(':')).toInt).toOption }

      // Because of Reasons, sometimes errors which aren't indicative of preemptions are treated as preemptions.
      val unsuccessfulStatusBuilder = errorCode match {
        case Status.ABORTED if jesCode.contains(GcpBatchAsyncBackendJobExecutionActor.JesPreemption) =>
          Preempted.apply _
        case Status.ABORTED
            if jesCode.contains(GcpBatchAsyncBackendJobExecutionActor.JesUnexpectedTermination) && wasPreemptible =>
          Preempted.apply _
        case Status.ABORTED if errorMessage.exists(_.contains(GcpBatchAsyncBackendJobExecutionActor.FailedV2Style)) =>
          Preempted.apply _
        case Status.UNKNOWN
            if errorMessage.exists(
              _.contains(GcpBatchAsyncBackendJobExecutionActor.FailedToStartDueToPreemptionSubstring)
            ) =>
          Preempted.apply _
        case Status.CANCELLED => Cancelled.apply _
        case _ => Failed.apply _
      }

      unsuccessfulStatusBuilder.apply(errorCode,
                                      jesCode,
                                      errorMessage.toList,
                                      eventList,
                                      machineType,
                                      zone,
                                      instanceName
      )
    }
  }

  case object DeletionInProgress extends RunStatus
  case object StateUnspecified extends RunStatus
  case object Unrecognized extends RunStatus

  final case class Failed(errorCode: Status,
                          jesCode: Option[Int],
                          errorMessages: List[String],
                          eventList: Seq[ExecutionEvent],
                          machineType: Option[String],
                          zone: Option[String],
                          instanceName: Option[String]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  /**
   * What Cromwell calls Aborted, PAPI calls Cancelled. This means the job was "cancelled" by the user
   */
  final case class Cancelled(errorCode: Status,
                             jesCode: Option[Int],
                             errorMessages: List[String],
                             eventList: Seq[ExecutionEvent],
                             machineType: Option[String],
                             zone: Option[String],
                             instanceName: Option[String]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Cancelled"
  }

  final case class Preempted(errorCode: Status,
                             jesCode: Option[Int],
                             errorMessages: List[String],
                             eventList: Seq[ExecutionEvent],
                             machineType: Option[String],
                             zone: Option[String],
                             instanceName: Option[String]
  ) extends UnsuccessfulRunStatus {
    override def toString = "Preempted"
  }
}
