package cromwell.backend.google.pipelines.common

import akka.actor.Props
import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.standard.pollmonitoring.{
  AsyncJobHasFinished,
  PollResultMessage,
  PollResultMonitorActor,
  ProcessThisPollResult
}
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

object PapiPollResultMonitorActor {
  def props(tellMetadataFn: Map[String, Any] => Unit,
            tellBardFn: (String, OffsetDateTime, OffsetDateTime, OffsetDateTime) => Unit
  ): Props = Props(new PapiPollResultMonitorActor(tellMetadataFn, tellBardFn))
}

class PapiPollResultMonitorActor(tellMetadataFn: Map[String, Any] => Unit,
                                 tellBardFn: (String, OffsetDateTime, OffsetDateTime, OffsetDateTime) => Unit
) extends PollResultMonitorActor[RunStatus] {

  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmStartTime => event.offsetDateTime
    }

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmEndTime => event.offsetDateTime
    }

  override def tellMetadata(metadata: Map[String, Any]): Unit = tellMetadataFn(metadata)
  override def tellBard(terminalStateName: String,
                        jobStart: OffsetDateTime,
                        vmStartTime: OffsetDateTime,
                        vmEndTime: OffsetDateTime
  ): Unit =
    tellBardFn(terminalStateName, jobStart, vmStartTime, vmEndTime)
  override def receive: Receive = {
    case message: PollResultMessage =>
      message match {
        case ProcessThisPollResult(pollResult: RunStatus) => processPollResult(pollResult)
        case ProcessThisPollResult(_) => println("Programmer error: Received Poll Result of unknown type.")
        case AsyncJobHasFinished(terminalStateName) => handleAsyncJobFinish(terminalStateName)
      }
    case _ =>
      println("Programmer error: Cost Helper received message of type other than CostPollingMessage")
  }
}
