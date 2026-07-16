package cromwell.backend.standard.pollmonitoring

import akka.actor.{Actor, ActorRef}
import cromwell.backend.validation._
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, Platform}
import cromwell.core.logging.JobLogger
import cromwell.services.cost.InstantiatedVmInfo
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime
trait PollResultMessage
case class ProcessThisPollResult[PollResultType](pollResult: PollResultType) extends PollResultMessage
case class AsyncJobHasFinished[PollResultType](pollResult: PollResultType) extends PollResultMessage

case class PollMonitorParameters(
  serviceRegistry: ActorRef,
  workflowDescriptor: BackendWorkflowDescriptor,
  jobDescriptor: BackendJobDescriptor,
  validatedRuntimeAttributes: ValidatedRuntimeAttributes,
  platform: Option[Platform],
  logger: JobLogger
)

/**
 * Processes poll results from backends and sends messages to other actors based on their contents.
 * Primarily concerned with reporting start times, end times, and cost data to the cromwell metadata service.
 */
trait PollResultMonitorActor[PollResultType] extends Actor {
  def params: PollMonitorParameters

  // Time that Cromwell (but not necessarily the cloud) started working on this job.
  def extractEarliestEventTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Time that the user VM started spending money.
  def extractStartTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Used to kick off a cost calculation
  def extractVmInfoFromRunState(pollStatus: PollResultType): Option[InstantiatedVmInfo]

  // Time that the user VM stopped spending money.
  def extractEndTimeFromRunState(pollStatus: PollResultType): Option[OffsetDateTime]

  // Function to emit metadata that is associated with a specific call attempt.
  def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    params.serviceRegistry.putMetadata(params.jobDescriptor.workflowDescriptor.id,
                                       Option(params.jobDescriptor.key),
                                       metadataKeyValues
    )
  }

  private var jobStartTime: Option[OffsetDateTime] =
    Option.empty
  private var vmStartTime: Option[OffsetDateTime] = Option.empty
  private var vmEndTime: Option[OffsetDateTime] = Option.empty
  protected var vmCostPerHour: Option[BigDecimal] = Option.empty

  def processPollResult(pollStatus: PollResultType): Unit = {
    // Make sure jobStartTime remains the earliest event time ever seen
    extractEarliestEventTimeFromRunState(pollStatus).foreach { earliestTime =>
      if (earliestTime.isBefore(jobStartTime.getOrElse(OffsetDateTime.now()))) {
        jobStartTime = Option(earliestTime)
      }
    }
    // If vm start time is reported and later than current start time, record it to metadata
    extractStartTimeFromRunState(pollStatus).foreach { start =>
      if (vmStartTime.isEmpty || start.isAfter(vmStartTime.get)) {
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
    // If we don't yet have a cost per hour and we can extract VM info, send a cost request to the catalog service.
    // We expect it to reply with an answer, which is handled in receive.
    // NB: Due to the nature of async code, we may send a few cost requests before we get a response back.
    if (vmCostPerHour.isEmpty) {
      extractVmInfoFromRunState(pollStatus).foreach(handleVmCostLookup)
    }
  }

  def handleVmCostLookup(vmInfo: InstantiatedVmInfo): Unit = ()

  def handleAsyncJobFinish(terminalStateName: String): Unit = ()
}
