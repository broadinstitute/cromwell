package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{ExecutionStatus, JobKey, JobOutputs, WorkflowId}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.services.MetadataServiceActor._
import cromwell.services._
import wdl4s._
import wdl4s.values.WdlValue

object EngineJobExecutionActor {

  /** States */
  sealed trait EngineJobExecutionActorState
  case object Pending extends EngineJobExecutionActorState
  case object CheckingJobStatus extends EngineJobExecutionActorState
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case class Start(jobKey: BackendJobDescriptorKey, data: WorkflowExecutionActorData, factory: BackendLifecycleActorFactory) extends EngineJobExecutionActorCommand
  case class Restart(jobKey: BackendJobDescriptorKey, data: WorkflowExecutionActorData) extends EngineJobExecutionActorCommand

  /** Responses */
  sealed trait EngineJobExecutionActorResponse
  case class JobSucceeded(jobKey: JobKey, returnCode: Option[Int], jobOutputs: JobOutputs)
  case class JobFailed(jobKey: JobKey, returnCode: Option[Int], reason: Throwable)

  object EngineJobExecutionActorData {
    def apply(currentActor: ActorRef) = new EngineJobExecutionActorData(Option(currentActor))
    def apply() = new EngineJobExecutionActorData(None)
  }
  case class EngineJobExecutionActorData(currentActor: Option[ActorRef])
}

class EngineJobExecutionActor(val workflowId: WorkflowId) extends LoggingFSM[EngineJobExecutionActorState, EngineJobExecutionActorData] with WorkflowLogging with ServiceRegistryClient {

  startWith(Pending, EngineJobExecutionActorData())

  when(Pending) {
    case Event(Start(jobKey, executionData, factory), _) =>
      val jobPreparationActorName = s"$workflowId-BackendPreparationActor-${jobKey.tag}"
      val jobPreparationActor = context.actorOf(JobPreparationActor.props(executionData, jobKey, factory), jobPreparationActorName)
      jobPreparationActor ! JobPreparationActor.Start
      goto(PreparingJob) using EngineJobExecutionActorData(jobPreparationActor)
  }

  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, actorProps), stateData) =>
      pushPreparedJobMetadata(jobDescriptor.key, jobDescriptor.inputs)
      val backendJobExecutionActor = context.actorOf(actorProps, buildJobExecutionActorName(jobDescriptor))
      backendJobExecutionActor ! ExecuteJobCommand
      pushRunningJobMetadata(jobDescriptor.key)
      goto(RunningJob) using EngineJobExecutionActorData(backendJobExecutionActor)
    case Event(response: BackendJobPreparationFailed, stateData) =>
      context.parent forward response
      context stop self
      stay()
  }

  when(RunningJob) {
    case Event(response: BackendJobExecutionResponse, stateData) =>
      context.parent forward response
      context stop self
      stay()
  }

  private def pushPreparedJobMetadata(jobKey: BackendJobDescriptorKey, jobInputs: Map[LocallyQualifiedName, WdlValue]) = {
    val inputEvents = jobInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKey(jobKey, s"${CallMetadataKeys.Inputs}")))
      case inputs =>
        inputs flatMap {
          case (inputName, inputValue) =>
            wdlValueToMetadataEvents(metadataKey(jobKey, s"${CallMetadataKeys.Inputs}:$inputName"), inputValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(inputEvents)
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }

  private def pushRunningJobMetadata(jobKey: JobKey) = {
    serviceRegistryActor ! PutMetadataAction(MetadataEvent(metadataKey(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Running)))
  }

  private def metadataKey(jobKey: JobKey, myKey: String) = MetadataKey(workflowId, Option(MetadataJobKey(jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)), myKey)
}
