package cromwell.backend.standard.costestimation

import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime



trait CostPollingHelper[A] {
  def extractStartTimeFromRunState(pollStatus: A): Option[OffsetDateTime]
  def extractEndTimeFromRunState(pollStatus: A): Option[OffsetDateTime]
  def calculateVmCostPerHour: BigDecimal
  def tellMetadata(metadata: Map[String, Any]): Unit

  var emittedStartTimeYet = false
  var emittedEndTimeYet = false
  def processPollResult(pollStatus: A): Unit = {
    if (!emittedStartTimeYet) {
      extractStartTimeFromRunState(pollStatus).foreach { start =>
        emittedStartTimeYet = true
        // NB: VM cost per hour will be emitted along with the start time.
        tellMetadata(Map(CallMetadataKeys.VmCostUsd -> calculateVmCostPerHour))
        tellMetadata(Map(CallMetadataKeys.VmStartTime -> start))
      }
    }
    if (!emittedEndTimeYet) {
      extractEndTimeFromRunState(pollStatus).foreach { end =>
        emittedEndTimeYet = true
        tellMetadata(Map(CallMetadataKeys.VmEndTime -> end))
      }
    }
  }
}

