package cromwell.backend.standard.pollmonitoring
import akka.actor.Actor
import cromwell.services.metadata.CallMetadataKeys
import java.time.OffsetDateTime
trait PollResultMessage
case class ProcessThisPollResult[PollResultType](pollResultType: PollResultType) extends PollResultMessage
case class AsyncJobHasFinished(terminalStateName: String) extends PollResultMessage

/**
 * Processes poll results from backends and sends messages to other actors based on their contents.
 * Primarily concerned with reporting start times, end times, and cost data to both the bard and cromwell metadata services.
 */
trait PollResultMonitorActor[PollResultType] extends Actor {
  // Time that Cromwell (but not necessarily the cloud) started working on this job.
  def extractEarliestEventTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Time that the user VM started spending money.
  def extractStartTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Time that the user VM stopped spending money.
  def extractEndTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Function to emit metadata that is associated with a specific call attempt.
  def tellMetadata(metadata: Map[String, Any]): Unit

  // Function that reports metrics to bard, called when a specific call attempt terminates.
  def tellBard(terminalStateName: String,
               jobStart: OffsetDateTime,
               vmStartTime: Option[OffsetDateTime],
               vmEndTime: OffsetDateTime
  ): Unit

  private var jobStartTime: Option[OffsetDateTime] =
    Option.empty
  private var vmStartTime: Option[OffsetDateTime] = Option.empty
  private var vmEndTime: Option[OffsetDateTime] = Option.empty

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
  def handleAsyncJobFinish(terminalStateName: String): Unit =
    jobStartTime.foreach(jobStart =>
      tellBard(terminalStateName, jobStart, vmStartTime, vmEndTime.getOrElse(OffsetDateTime.now()))
    )
}
