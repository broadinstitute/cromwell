package cromwell.backend.standard.costestimation
import akka.actor.Actor
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

// Emits metadata based on backend poll results. Limits the amount of metadata emitted.
trait CostPollingHelper[A] extends Actor {

  // Should be overridden to return the time that the user VM starts spending money.
  // None if it can't be ascertained from the provided status.
  def extractStartTimeFromRunState(pollStatus: A): Option[OffsetDateTime]

  // Should be overridden to return the time that the user VM stops spending money.
  // None if it can't be ascertained from the provided status.
  def extractEndTimeFromRunState(pollStatus: A): Option[OffsetDateTime]
  // Function to emit metadata that is associated with a specific call attempt.

  def extractVmCostPerHourFromRunState(pollStatus: A): Option[BigDecimal]
  def tellMetadata(metadata: Map[String, Any]): Unit

  var vmStartTime: Option[OffsetDateTime] = Option.empty
  var vmEndTime: Option[OffsetDateTime] = Option.empty
  var vmCostPerHour: Option[BigDecimal] = Option.empty
  def processPollResult(pollStatus: A): Unit = {
    if (vmStartTime.isEmpty) {
      extractStartTimeFromRunState(pollStatus).foreach { start =>
        vmStartTime = Some(start)
        tellMetadata(Map(CallMetadataKeys.VmStartTime -> start))
      }
    }

    extractEndTimeFromRunState(pollStatus).foreach { end =>
      if (vmEndTime.isEmpty || end.isAfter(vmEndTime.get)) {
        vmEndTime = Some(end)
        tellMetadata(Map(CallMetadataKeys.VmEndTime -> end))
      }
    }

    if (vmCostPerHour.isEmpty) {
      extractVmCostPerHourFromRunState(pollStatus).foreach { costPerHour =>
        vmCostPerHour = Option(costPerHour)
        tellMetadata(Map(CallMetadataKeys.VmCostPerHour -> costPerHour))
      }
    }
  }
}
