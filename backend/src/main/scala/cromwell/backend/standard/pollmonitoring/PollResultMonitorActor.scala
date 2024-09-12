package cromwell.backend.standard.pollmonitoring
import akka.actor.Actor
import cromwell.services.metadata.CallMetadataKeys
import java.time.OffsetDateTime
trait PollResultMessage
case class ProcessThisPollResult[PollResultType](pollResultType: PollResultType) extends PollResultMessage
case class AsyncJobHasFinished(terminalStateName: String) extends PollResultMessage

// Processes poll results from backends and sends messages to other actors based on their contents.
// Primarily concerned with reporting start times, end times, and cost data to both the bard service and cromwell metadata.
trait PollResultMonitorActor[PollResultType] extends Actor {

  // Should be overridden to return the time that the user VM starts spending money.
  // None if it can't be ascertained from the provided status.
  def extractStartTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Should be overridden to return the time that the user VM stops spending money.
  // None if it can't be ascertained from the provided status.
  def extractEndTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Function to emit metadata that is associated with a specific call attempt.
  def tellMetadata(metadata: Map[String, Any]): Unit

  // Function that reports metrics to bard, called when a specific call attempt terminates.
  def tellBard(terminalStateName: String,
               jobStart: OffsetDateTime,
               vmStartTime: OffsetDateTime,
               vmEndTime: OffsetDateTime
  ): Unit

  var vmStartTime: Option[OffsetDateTime] = Option.empty
  var vmEndTime: Option[OffsetDateTime] = Option.empty

  def processPollResult(pollStatus: PollResultType): Unit = {
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
  }

  def handleAsyncJobFinish(terminalStateName: String): Unit =
    tellBard(terminalStateName, vmStartTime.get, vmStartTime.get, vmEndTime.get)
}
