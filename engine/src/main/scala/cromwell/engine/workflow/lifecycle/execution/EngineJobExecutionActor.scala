package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.WorkflowId
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.jobstore._

object EngineJobExecutionActor {
  /** States */
  sealed trait EngineJobExecutionActorState
  case object Pending extends EngineJobExecutionActorState
  case object CheckingJobStatus extends EngineJobExecutionActorState
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case class Start(jobKey: BackendJobDescriptorKey) extends EngineJobExecutionActorCommand
  case class Restart(jobKey: BackendJobDescriptorKey) extends EngineJobExecutionActorCommand

  case class JobRunning(jobDescriptor: BackendJobDescriptor, backendJobExecutionActor: ActorRef)

  protected object EngineJobExecutionActorData {
    def apply(currentActor: ActorRef) = new EngineJobExecutionActorData(Option(currentActor))
    def apply() = new EngineJobExecutionActorData(None)
  }

  protected case class EngineJobExecutionActorData(currentActor: Option[ActorRef])

  def props(workflowId: WorkflowId, executionData: WorkflowExecutionActorData,
            factory: BackendLifecycleActorFactory, initializationData: Option[BackendInitializationData]) = {
    Props(new EngineJobExecutionActor(workflowId, executionData, factory, initializationData))
  }
}

class EngineJobExecutionActor(val workflowId: WorkflowId,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData]) extends LoggingFSM[EngineJobExecutionActorState, EngineJobExecutionActorData] with WorkflowLogging {

  startWith(Pending, EngineJobExecutionActorData())

  when(Pending) {
    case Event(Start(jobKey), _) =>
      val jobPreparationActorName = s"$workflowId-BackendPreparationActor-${jobKey.tag}"
      val jobPreparationActor = context.actorOf(JobPreparationActor.props(executionData, jobKey, factory, initializationData), jobPreparationActorName)
      jobPreparationActor ! JobPreparationActor.Start
      goto(PreparingJob) using EngineJobExecutionActorData(jobPreparationActor)
    case Event(Restart(jobKey), _) =>
      context.actorOf(JobStoreReader.props()) ! QueryJobCompletion(jobKey)
      goto(CheckingJobStatus)
  }

  when(CheckingJobStatus) {
    case Event(JobNotComplete(jobKey), _) =>
      self ! Start(jobKey)
      goto(Pending)
    case Event(JobComplete(jobKey, jobResult), _) =>
      jobResult match {
        case JobResultSuccess(returnCode, jobOutputs) =>
          context.parent ! new SucceededResponse(jobKey, returnCode, jobOutputs)
          context stop self
          stay()
        case JobResultFailure(returnCode, reason) =>
          context.parent ! new FailedNonRetryableResponse(jobKey, reason, returnCode)
          context stop self
          stay()
      }
  }

  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, actorProps), stateData) =>
      val backendJobExecutionActor = context.actorOf(actorProps, buildJobExecutionActorName(jobDescriptor))
      backendJobExecutionActor ! ExecuteJobCommand
      context.parent ! JobRunning(jobDescriptor, backendJobExecutionActor)
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
  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }
}
