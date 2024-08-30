package cromwell.backend.google.pipelines.v2beta.api

import cromwell.backend.standard.costestimation.CostPollingHelper
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.google.pipelines.common.api.RunStatus.{Running, TerminalRunStatus}
import java.time.OffsetDateTime

class PapiCostPollingHelper(tellMetadataFn: Map[String, Any] => Unit) extends CostPollingHelper[RunStatus] {

  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] = pollStatus match {
    case runningStatus: Running => runningStatus.vmStartTime
    case _ => Option.empty
  }

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] = pollStatus match {
    case terminalStatus: TerminalRunStatus => terminalStatus.vmEndTime
    case _ => Option.empty
  }

  override def calculateVmCostPerHour: BigDecimal = 3.50 // TODO

  override def tellMetadata(metadata: Map[String, Any]): Unit = tellMetadataFn(metadata)
}
