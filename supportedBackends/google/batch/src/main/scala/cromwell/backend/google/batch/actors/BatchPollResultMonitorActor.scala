package cromwell.backend.google.batch.actors

import akka.actor.{ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, Platform}
import cromwell.backend.google.batch.models.RunStatus
import cromwell.backend.standard.pollmonitoring.{
  AsyncJobHasFinished,
  PollMonitorParameters,
  PollResultMessage,
  PollResultMonitorActor,
  ProcessThisPollResult
}
import cromwell.backend.validation.ValidatedRuntimeAttributes
import cromwell.core.logging.JobLogger
import cromwell.services.cost.{GcpCostLookupRequest, GcpCostLookupResponse, InstantiatedVmInfo}
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

object BatchPollResultMonitorActor {
  def props(serviceRegistry: ActorRef,
            workflowDescriptor: BackendWorkflowDescriptor,
            jobDescriptor: BackendJobDescriptor,
            runtimeAttributes: ValidatedRuntimeAttributes,
            platform: Option[Platform],
            logger: JobLogger
  ): Props = Props(
    new BatchPollResultMonitorActor(
      PollMonitorParameters(serviceRegistry, workflowDescriptor, jobDescriptor, runtimeAttributes, platform, logger)
    )
  )
}

class BatchPollResultMonitorActor(pollMonitorParameters: PollMonitorParameters)
    extends PollResultMonitorActor[RunStatus] {

  override def extractEarliestEventTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.minByOption(_.offsetDateTime).map(e => e.offsetDateTime)

  // gets the newest VmStartTime/VmEndTime because the user is only charged for the last time a job transitions to SCHEDULED state
  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList
      .filter(_.name == CallMetadataKeys.VmStartTime)
      .maxByOption(_.offsetDateTime)
      .map(e => e.offsetDateTime)

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList
      .filter(_.name == CallMetadataKeys.VmEndTime)
      .maxByOption(_.offsetDateTime)
      .map(e => e.offsetDateTime)

  override def extractVmInfoFromRunState(pollStatus: RunStatus): Option[InstantiatedVmInfo] =
    pollStatus.instantiatedVmInfo

  override def handleVmCostLookup(vmInfo: InstantiatedVmInfo) = {
    val request = GcpCostLookupRequest(vmInfo, self)
    params.serviceRegistry ! request
  }

  def handleCostResponse(costLookupResponse: GcpCostLookupResponse): Unit =
    if (vmCostPerHour.isEmpty) { // Optimization to avoid processing responses after we've received a valid one.
      val cost = costLookupResponse.calculatedCost match {
        case Valid(c) =>
          params.logger.info(s"vmCostPerHour for ${costLookupResponse.vmInfo} is ${c}")
          c
        case Invalid(errors) =>
          params.logger.error(
            s"Failed to calculate VM cost per hour for ${costLookupResponse.vmInfo}. ${errors.toList.mkString(", ")}"
          )
          BigDecimal(-1)
      }
      vmCostPerHour = Option(cost)
      tellMetadata(Map(CallMetadataKeys.VmCostPerHour -> cost))
    }

  override def receive: Receive = {
    case costResponse: GcpCostLookupResponse => handleCostResponse(costResponse)
    case message: PollResultMessage =>
      message match {
        case ProcessThisPollResult(pollResult: RunStatus) => processPollResult(pollResult)
        case ProcessThisPollResult(result) =>
          params.logger.error(
            s"Programmer error: Received Poll Result of unknown type. Expected ${RunStatus.getClass.getSimpleName} but got ${result.getClass.getSimpleName}."
          )

        case AsyncJobHasFinished(pollResult: RunStatus) => handleAsyncJobFinish(pollResult.getClass.getSimpleName)
        case AsyncJobHasFinished(result) =>
          params.logger.error(
            s"Programmer error: Received Poll Result of unknown type. Expected ${AsyncJobHasFinished.getClass.getSimpleName} but got ${result.getClass.getSimpleName}."
          )

      }
    case _ =>
      params.logger.error(
        s"Programmer error: Cost Helper received message of type other than CostPollingMessage"
      )

  }

  override def params: PollMonitorParameters = pollMonitorParameters

}
