package cromwell.backend.google.batch.models

//import com.google.cloud.batch.v1.JobStatus
import cromwell.backend.google.batch.actors.GcpBatchAsyncBackendJobExecutionActor
import cromwell.core.ExecutionEvent
import io.grpc.Status
//import org.slf4j.{Logger, LoggerFactory}

import scala.util.Try

sealed trait RunStatus

object RunStatus {

//  private val log: Logger = LoggerFactory.getLogger(RunStatus.toString)
//  def fromJobStatus(status: JobStatus.State): RunStatus = {
//    log.error(s"REmove me: $status")
//    ???
//  }

  // TODO: Fix this

//  def fromJobStatus(status: JobStatus.State): RunStatus = status match {
//    case JobStatus.State.QUEUED =>
//      log.info("job queued")
//      Running
//    case JobStatus.State.SCHEDULED =>
//      log.info("job scheduled")
//      Running
//    case JobStatus.State.RUNNING =>
//      log.info("job running")
//      Running
//    case JobStatus.State.SUCCEEDED =>
//      log.info("job scheduled")
//      Succeeded(List(ExecutionEvent("complete in GCP Batch"))) // update to more specific
//    case JobStatus.State.FAILED =>
//      log.info("job failed")
//      Failed(List.empty)
//    case JobStatus.State.DELETION_IN_PROGRESS =>
//      log.info("deletion in progress")
//      DeletionInProgress
//    case JobStatus.State.STATE_UNSPECIFIED =>
//      log.info("state unspecified")
//      StateUnspecified
//    case JobStatus.State.UNRECOGNIZED =>
//      log.info("state unrecognized")
//      Unrecognized
//    case _ =>
//      log.info(s"job status not matched: $status")
//      Running
//  }

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

//  case class Succeeded(override val eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
//    override def toString = "Succeeded"
//  }

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
