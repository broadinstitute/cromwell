package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
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

  def props(executionData: WorkflowExecutionActorData, factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData]) = {
    Props(new EngineJobExecutionActor(executionData, factory, initializationData))
  }

  /** Data */
  case class EngineJobExecutionActorData(restarting: Boolean = false)
}

class EngineJobExecutionActor(executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData]) extends LoggingFSM[EngineJobExecutionActorState, EngineJobExecutionActorData] with WorkflowLogging {

  override val workflowId = executionData.workflowDescriptor.id

  startWith(Pending, EngineJobExecutionActorData())

  when(Pending) {
    case Event(Start(jobKey), _) =>
      prepareJob(jobKey) using EngineJobExecutionActorData(restarting = false)
    case Event(Restart(jobKey), _) =>
      context.actorOf(JobStoreReader.props()) ! QueryJobCompletion(jobKey)
      goto(CheckingJobStatus) using EngineJobExecutionActorData(restarting = true)
  }

  when(CheckingJobStatus) {
    case Event(JobNotComplete(jobKey), _) =>
      prepareJob(jobKey)
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
      if (stateData.restarting) backendJobExecutionActor ! RecoverJobCommand
      else backendJobExecutionActor ! ExecuteJobCommand
      context.parent ! JobRunning(jobDescriptor, backendJobExecutionActor)
      goto(RunningJob)
    case Event(response: BackendJobPreparationFailed, _) =>
      context.parent forward response
      context stop self
      stay()
  }

  when(RunningJob) {
    case Event(response: BackendJobExecutionResponse, _) =>
      context.parent forward response
      context stop self
      stay()
  }

  def prepareJob(jobKey: BackendJobDescriptorKey) = {
    val jobPreparationActorName = s"$workflowId-BackendPreparationActor-${jobKey.tag}"
    val jobPreparationActor = context.actorOf(JobPreparationActor.props(executionData, jobKey, factory, initializationData), jobPreparationActorName)
    jobPreparationActor ! JobPreparationActor.Start
    goto(PreparingJob)
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }
}
