package cromwell.backend.standard.pollmonitoring
import akka.actor.{Actor, ActorRef}
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, Platform}
import cromwell.backend.validation.{
  CpuValidation,
  DockerValidation,
  MemoryValidation,
  RuntimeAttributesValidation,
  ValidatedRuntimeAttributes
}
import cromwell.core.logging.JobLogger
import cromwell.services.cost.{GcpCostLookupRequest, GcpCostLookupResponse, InstantiatedVmInfo}
import cromwell.services.metadata.CallMetadataKeys
import cromwell.services.metrics.bard.BardEventing.BardEventRequest
import cromwell.services.metrics.bard.model.TaskSummaryEvent
import wdl4s.parser.MemoryUnit

import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
trait PollResultMessage
case class ProcessThisPollResult[PollResultType](pollResult: PollResultType) extends PollResultMessage
case class AsyncJobHasFinished[PollResultType](pollResult: PollResultType) extends PollResultMessage

case class PollMonitorParameters(
  serviceRegistry: ActorRef,
  workflowDescriptor: BackendWorkflowDescriptor,
  jobDescriptor: BackendJobDescriptor,
  validatedRuntimeAttributes: ValidatedRuntimeAttributes,
  platform: Option[Platform],
  logger: Option[JobLogger]
)

/**
 * Processes poll results from backends and sends messages to other actors based on their contents.
 * Primarily concerned with reporting start times, end times, and cost data to both the bard and cromwell metadata services.
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

  // Function that reports metrics to bard, called when a specific call attempt terminates.
  def tellBard(terminalStateName: String,
               jobStart: OffsetDateTime,
               vmStartTime: Option[OffsetDateTime],
               vmEndTime: OffsetDateTime
  ): Unit = {
    val validatedRuntimeAttributes = params.validatedRuntimeAttributes
    val serviceRegistryActor = params.serviceRegistry
    val workflowDescriptor = params.workflowDescriptor
    val jobDescriptor = params.jobDescriptor
    val platform = params.platform.map(_.runtimeKey)
    val dockerImage =
      RuntimeAttributesValidation.extractOption(DockerValidation.instance, validatedRuntimeAttributes)
    val cpus = RuntimeAttributesValidation.extract(CpuValidation.instance, validatedRuntimeAttributes).value
    val memory = RuntimeAttributesValidation
      .extract(MemoryValidation.instance(), validatedRuntimeAttributes)
      .to(MemoryUnit.Bytes)
      .amount
    serviceRegistryActor ! BardEventRequest(
      TaskSummaryEvent(
        workflowDescriptor.id.id,
        workflowDescriptor.possibleParentWorkflowId.map(_.id),
        workflowDescriptor.rootWorkflowId.id,
        jobDescriptor.key.tag,
        jobDescriptor.key.call.fullyQualifiedName,
        jobDescriptor.key.index,
        jobDescriptor.key.attempt,
        terminalStateName,
        platform,
        dockerImage,
        cpus,
        memory,
        jobStart.toString,
        vmStartTime.map(startTime => startTime.toString),
        vmEndTime.toString,
        jobStart.until(vmEndTime, ChronoUnit.SECONDS),
        vmStartTime.map(start => start.until(vmEndTime, ChronoUnit.SECONDS))
      )
    )
  }

  private var jobStartTime: Option[OffsetDateTime] =
    Option.empty
  private var vmStartTime: Option[OffsetDateTime] = Option.empty
  private var vmEndTime: Option[OffsetDateTime] = Option.empty
  private var vmCostPerHour: Option[BigDecimal] = Option.empty

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
    // If we don't yet have a cost per hour and we can extract VM info, send a cost request to the catalog service.
    // We expect it to reply with an answer, which is handled in receive.
    // NB: Due to the nature of async code, we may send a few cost requests before we get a response back.
    if (vmCostPerHour.isEmpty) {
      val instantiatedVmInfo = extractVmInfoFromRunState(pollStatus)
      instantiatedVmInfo.foreach { vmInfo =>
        println(s"Requesting cost info for: ${vmInfo}")
        val request = GcpCostLookupRequest(vmInfo, self)
        params.serviceRegistry ! request
      }
    }
  }

  // When a job finishes, the bard actor needs to know about the timing in order to record metrics.
  // Cost related metadata should already have been handled in processPollResult.
  def handleAsyncJobFinish(terminalStateName: String): Unit =
    jobStartTime.foreach(jobStart =>
      tellBard(
        terminalStateName = terminalStateName,
        jobStart = jobStart,
        vmStartTime = vmStartTime,
        vmEndTime = vmEndTime.getOrElse(OffsetDateTime.now())
      )
    )

  def handleCostResponse(costLookupResponse: GcpCostLookupResponse): Unit = {
    println(s"Handling Cost Response from Catalog Service: ${costLookupResponse}")
    if (vmCostPerHour.isEmpty) { // Optimization to avoid processing responses after we've received a valid one.
      val cost = costLookupResponse.calculatedCost match {
        case Valid(c) => c
        case Invalid(errors) =>
          // TODO contextualizeErrors
          params.logger.foreach(_.error(s"Failed to calculate VM cost per hour. ${errors.toList.mkString(", ")}"))
          BigDecimal(-1)
      }
      vmCostPerHour = Option(cost)
      tellMetadata(Map(CallMetadataKeys.VmCostPerHour -> cost))
    }
  }
}
