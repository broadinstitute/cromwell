package cromwell.backend.google.pipelines.batch

import cromwell.core.ExecutionEvent

sealed trait GcpBatchRunStatus

object GcpBatchRunStatus {
    //case object Initializing extends GcpBatchRunStatus
    //case object AwaitingCloudQuota extends RunStatus
    case object Running extends GcpBatchRunStatus {
      def isTerminal = false
    }
    case object Complete extends GcpBatchRunStatus {
      def isTerminal = true
    }

    sealed trait TerminalRunStatus extends GcpBatchRunStatus {
      def eventList: Seq[ExecutionEvent]
      def machineType: Option[String]
      def zone: Option[String]
      def instanceName: Option[String]


    }

    case class Running()

    case class Complete()

    case class Success(eventList: Seq[ExecutionEvent])
                       //machineType: Option[String],
                       //zone: Option[String])
                       //instanceName: Option[String]) extends TerminalRunStatus {
        override def toString = "Success"

}
