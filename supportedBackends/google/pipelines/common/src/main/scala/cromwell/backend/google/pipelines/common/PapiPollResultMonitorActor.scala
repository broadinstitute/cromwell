package cromwell.backend.google.pipelines.common

import akka.actor.{ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.standard.pollmonitoring._
import cromwell.backend.validation.ValidatedRuntimeAttributes
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, Platform}
import cromwell.core.logging.JobLogger
import cromwell.services.cost.{GcpCostLookupRequest, GcpCostLookupResponse, InstantiatedVmInfo}
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

object PapiPollResultMonitorActor {
  def props(serviceRegistry: ActorRef,
            workflowDescriptor: BackendWorkflowDescriptor,
            jobDescriptor: BackendJobDescriptor,
            runtimeAttributes: ValidatedRuntimeAttributes,
            platform: Option[Platform],
            logger: JobLogger
  ): Props = Props(
    new PapiPollResultMonitorActor(
      PollMonitorParameters(serviceRegistry, workflowDescriptor, jobDescriptor, runtimeAttributes, platform, logger)
    )
  )
}

class PapiPollResultMonitorActor(parameters: PollMonitorParameters) extends PollResultMonitorActor[RunStatus] {

  override def extractEarliestEventTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.minByOption(_.offsetDateTime).map(e => e.offsetDateTime)

  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmStartTime => event.offsetDateTime
    }

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmEndTime => event.offsetDateTime
    }

  override def extractVmInfoFromRunState(pollStatus: RunStatus): Option[InstantiatedVmInfo] =
    pollStatus.instantiatedVmInfo

  override def handleVmCostLookup(vmInfo: InstantiatedVmInfo) = {
    val request = GcpCostLookupRequest(vmInfo, self)
    params.serviceRegistry ! request
  }

  override def params: PollMonitorParameters = parameters

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
    case unexpected =>
      params.logger.error(
        s"Programmer error: Cost Helper received message of unexpected type. Was ${unexpected.getClass.getSimpleName}."
      )

  }
}
