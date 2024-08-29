package cromwell.backend.google.pipelines.common.api

import _root_.io.grpc.Status
import cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor
import cromwell.core.ExecutionEvent

import java.time.OffsetDateTime
import scala.util.Try

sealed trait RunStatus

object RunStatus {
  case object Initializing extends RunStatus
  case object AwaitingCloudQuota extends RunStatus
  case class Running(vmStartTime: Option[OffsetDateTime]) extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
    def machineType: Option[String]
    def zone: Option[String]
    def instanceName: Option[String]
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

  case class Success(eventList: Seq[ExecutionEvent],
                     machineType: Option[String],
                     zone: Option[String],
                     instanceName: Option[String]
  ) extends TerminalRunStatus {
    override def toString = "Success"
  }

  object UnsuccessfulRunStatus {

    // TODO: Dead code alert. Functional code only ever calls this with status `UNKNOWN`.
    //
    // Seems to have been replaced with:
    //   - cromwell.backend.google.pipelines.v2beta.api.request.ErrorReporter#toUnsuccessfulRunStatus
    //   - cromwell.backend.google.batch.models.RunStatus#fromJobStatus
    // There are useful tests in `PipelinesApiAsyncBackendJobExecutionActorSpec`
    // that test other things and happen to rely on this method, so eventually
    // delete it with the rest of Life Sciences. GCP Batch does not use the
    // `PipelinesApiAsyncBackendJobExecutionActor` at all.
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
        case Status.ABORTED if jesCode.contains(PipelinesApiAsyncBackendJobExecutionActor.JesPreemption) =>
          Preempted.apply _
        case Status.ABORTED
            if jesCode.contains(PipelinesApiAsyncBackendJobExecutionActor.JesUnexpectedTermination) && wasPreemptible =>
          Preempted.apply _
        case Status.ABORTED
            if errorMessage.exists(_.contains(PipelinesApiAsyncBackendJobExecutionActor.FailedV2Style)) =>
          Preempted.apply _
        case Status.UNKNOWN
            if errorMessage.exists(
              _.contains(PipelinesApiAsyncBackendJobExecutionActor.FailedToStartDueToPreemptionSubstring)
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

  /**
    * This should NOT happen, but occasionally we see Life Sciences fail jobs with
    * as FAILED_PRECONDITION and a message that contains "no available zones" or similar. (WX-1625)
    */
  final case class QuotaFailed(errorCode: Status,
                               jesCode: Option[Int],
                               errorMessages: List[String],
                               eventList: Seq[ExecutionEvent],
                               machineType: Option[String],
                               zone: Option[String],
                               instanceName: Option[String]
  ) extends UnsuccessfulRunStatus {
    override def toString = "QuotaFailed"
  }
}
