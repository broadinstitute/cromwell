package cromwell.backend.standard.costestimation
import cromwell.services.metadata.CallMetadataKeys
import java.time.OffsetDateTime

trait CostPollingHelper[A] {
  def extractStartTimeFromRunState(pollStatus: A): Option[OffsetDateTime]
  def extractEndTimeFromRunState(pollStatus: A): Option[OffsetDateTime]
  def calculateVmCostPerHour: BigDecimal
  def tellMetadata(metadata: Map[String, Any]): Unit

  var vmStartTime: Option[OffsetDateTime] = Option.empty
  var vmEndTime: Option[OffsetDateTime] = Option.empty
  def processPollResult(pollStatus: A): Unit = {
    if (vmStartTime.isEmpty) {
      extractStartTimeFromRunState(pollStatus).foreach { start =>
        vmStartTime = Some(start)
        // NB: VM cost per hour will be emitted along with the start time.
        tellMetadata(Map(CallMetadataKeys.VmCostUsd -> calculateVmCostPerHour))
        tellMetadata(Map(CallMetadataKeys.VmStartTime -> start))
      }
    }
    if (vmEndTime.isEmpty) {
      extractEndTimeFromRunState(pollStatus).foreach { end =>
        vmEndTime = Some(end)
        tellMetadata(Map(CallMetadataKeys.VmEndTime -> end))
      }
    }
  }
}
