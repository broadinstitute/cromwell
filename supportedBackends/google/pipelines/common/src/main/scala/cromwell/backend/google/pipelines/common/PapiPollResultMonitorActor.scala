package cromwell.backend.google.pipelines.common

import akka.actor.{ActorRef, Props}
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.standard.pollmonitoring._
import cromwell.backend.validation.ValidatedRuntimeAttributes
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, Platform}
import cromwell.core.logging.JobLogger
import cromwell.services.cost.{GcpCostLookupResponse, InstantiatedVmInfo}
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
      PollMonitorParameters(serviceRegistry,
                            workflowDescriptor,
                            jobDescriptor,
                            runtimeAttributes,
                            platform,
                            Option(logger)
      )
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

  override def params: PollMonitorParameters = parameters

  override def receive: Receive = {
    case costResponse: GcpCostLookupResponse => handleCostResponse(costResponse)
    case message: PollResultMessage =>
      message match {
        case ProcessThisPollResult(pollResult: RunStatus) => processPollResult(pollResult)
        case ProcessThisPollResult(result) =>
          params.logger.foreach(logger =>
            logger.error(
              s"Programmer error: Received Poll Result of unknown type. Expected ${RunStatus.getClass.getSimpleName} but got ${result.getClass.getSimpleName}."
            )
          )
        case AsyncJobHasFinished(pollResult: RunStatus) => handleAsyncJobFinish(pollResult.getClass.getSimpleName)
        case AsyncJobHasFinished(result) =>
          params.logger.foreach(logger =>
            logger.error(
              s"Programmer error: Received Poll Result of unknown type. Expected ${AsyncJobHasFinished.getClass.getSimpleName} but got ${result.getClass.getSimpleName}."
            )
          )
      }
    case _ =>
      params.logger.foreach(logger =>
        logger.error(
          s"Programmer error: Cost Helper received message of type other than CostPollingMessage"
        )
      )
  }
}
