package cromwell.backend.google.pipelines.batch

import cromwell.core.ExecutionEvent

sealed trait GcpBatchRunStatus

object GcpBatchRunStatus {
    case object Initializing extends GcpBatchRunStatus
    //case object AwaitingCloudQuota extends RunStatus
    case object Running extends GcpBatchRunStatus

    case class Success(eventList: Seq[ExecutionEvent],
                       machineType: Option[String],
                       zone: Option[String])
                       //instanceName: Option[String]) extends TerminalRunStatus {
        override def toString = "Success"

}
