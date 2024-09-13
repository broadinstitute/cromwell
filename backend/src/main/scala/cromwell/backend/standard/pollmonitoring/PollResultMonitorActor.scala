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

  // Should be overridden to return the earliest time that anything has happened. Used to determine when Cromwell (but not necessarily the cloud)
  // started working on this job.
  def extractEarliestEventTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

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
               vmStartTime: Option[OffsetDateTime],
               vmEndTime: OffsetDateTime
  ): Unit

  var jobStartTime: Option[OffsetDateTime] =
    Option.empty // earliest recorded time that this job has done anything at all
  var vmStartTime: Option[OffsetDateTime] = Option.empty // time that the VM starts spending money
  var vmEndTime: Option[OffsetDateTime] = Option.empty // time that the VM stopped spending money, falling back to now.

  def processPollResult(pollStatus: PollResultType): Unit = {
    // Make sure jobStartTime remains the earliest event time ever seen
    extractEarliestEventTimeFromRunState(pollStatus).foreach { earliestTime =>
      if (earliestTime.isBefore(jobStartTime.getOrElse(OffsetDateTime.now()))) {
        jobStartTime = Option(earliestTime)
      }
    }
    // If vm start time is reported, record it to metadata and stop trying
    if (vmStartTime.isEmpty) {
      extractStartTimeFromRunState(pollStatus).foreach { start =>
        vmStartTime = Option(start)
        tellMetadata(Map(CallMetadataKeys.VmStartTime -> start))
      }
    }
    // If vm end time is reported, (or for some weird reason we see an end time after our recorded one),
    // record it to metadata.
    extractEndTimeFromRunState(pollStatus).foreach { end =>
      if (vmEndTime.isEmpty || end.isAfter(vmEndTime.get)) {
        vmEndTime = Option(end)
        tellMetadata(Map(CallMetadataKeys.VmEndTime -> end))
      }
    }
  }

  // When a job finishes, the bard actor needs to know about the timing in order to record metrics.
  // Cost related metadata should already have been handled in processPollResult.
  def handleAsyncJobFinish(terminalStateName: String): Unit = {
    jobStartTime.foreach(jobStart => tellBard(terminalStateName, jobStart, vmStartTime, vmEndTime.getOrElse(OffsetDateTime.now())))
  }
}
