package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.standard.costestimation.CostPollingHelper
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

class PapiCostPollingHelper(tellMetadataFn: Map[String, Any] => Unit) extends CostPollingHelper[RunStatus] {

  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmStartTime => event.offsetDateTime
    }

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmEndTime => event.offsetDateTime
    }

  override def calculateVmCostPerHour: BigDecimal = 3.50 // TODO

  override def tellMetadata(metadata: Map[String, Any]): Unit = tellMetadataFn(metadata)
}
