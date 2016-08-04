package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.logging.WorkflowLogging
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.{BackendJobPreparationFailed, BackendJobPreparationSucceeded}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.{Pending => _, _}

object EngineJobExecutionActor {
  /** States */
  sealed trait EngineJobExecutionActorState
  case object Pending extends EngineJobExecutionActorState
  case object CheckingJobStatus extends EngineJobExecutionActorState
  case object PreparingJob extends EngineJobExecutionActorState
  case object RunningJob extends EngineJobExecutionActorState
  case object AwaitingJobStoreAck extends EngineJobExecutionActorState

  /** Commands */
  sealed trait EngineJobExecutionActorCommand
  case object Execute extends EngineJobExecutionActorCommand

  final case class JobRunning(jobDescriptor: BackendJobDescriptor, backendJobExecutionActor: ActorRef)

  def props(jobDescriptorKey: BackendJobDescriptorKey,
            executionData: WorkflowExecutionActorData,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            restarting: Boolean,
            serviceRegistryActor: ActorRef,
            jobStoreActor: ActorRef) = {
    Props(new EngineJobExecutionActor(jobDescriptorKey,
      executionData,
      factory,
      initializationData,
      restarting,
      serviceRegistryActor,
      jobStoreActor)).withDispatcher(EngineDispatcher)
  }
}

class EngineJobExecutionActor(jobKey: BackendJobDescriptorKey,
                              executionData: WorkflowExecutionActorData,
                              factory: BackendLifecycleActorFactory,
                              initializationData: Option[BackendInitializationData],
                              restarting: Boolean,
                              serviceRegistryActor: ActorRef,
                              jobStoreActor: ActorRef) extends LoggingFSM[EngineJobExecutionActorState, Unit] with WorkflowLogging {

  override val workflowId = executionData.workflowDescriptor.id

  startWith(Pending, ())

  when(Pending) {
    case Event(Execute, _) =>
      if (restarting) {
        val jobStoreKey = jobKey.toJobStoreKey(workflowId)
        jobStoreActor ! QueryJobCompletion(jobStoreKey)
        goto(CheckingJobStatus)
      } else {
        prepareJob(jobKey)
      }
  }

  when(CheckingJobStatus) {
    case Event(JobNotComplete, _) =>
      prepareJob(jobKey)
    case Event(JobComplete(jobResult), _) =>
      jobResult match {
        case JobResultSuccess(returnCode, jobOutputs) =>
          context.parent ! SucceededResponse(jobKey, returnCode, jobOutputs)
          context stop self
          stay()
        case JobResultFailure(returnCode, reason, false) =>
          context.parent ! FailedNonRetryableResponse(jobKey, reason, returnCode)
          context stop self
          stay()
        case JobResultFailure(returnCode, reason, true) =>
          context.parent ! FailedRetryableResponse(jobKey, reason, returnCode)
          context stop self
          stay()
      }
    case Event(f: JobStoreReadFailure, _) =>
      log.error(f.reason, "Error reading from JobStore for " + jobKey)
      // Escalate
      throw new RuntimeException(f.reason)
  }

  when(PreparingJob) {
    case Event(BackendJobPreparationSucceeded(jobDescriptor, actorProps), stateData) =>
      val backendJobExecutionActor = context.actorOf(actorProps, buildJobExecutionActorName(jobDescriptor))
      val message = if (restarting) RecoverJobCommand else ExecuteJobCommand
      backendJobExecutionActor ! message
      context.parent ! JobRunning(jobDescriptor, backendJobExecutionActor)
      goto(RunningJob)
    case Event(response: BackendJobPreparationFailed, _) =>
      context.parent forward response
      context stop self
      stay()
  }

  when(RunningJob) {
    case Event(response: BackendJobExecutionResponse, _) =>
      saveJobCompletionResult(response)
      context.parent forward response
      goto(AwaitingJobStoreAck)
  }

  when(AwaitingJobStoreAck) {
    case Event(JobStoreWriteSuccess(_), _) =>
      context stop self
      stay()
    case Event(JobStoreWriteFailure(_, t), _) =>
      // This is moderately bad: If we can't write this result to the database, and then undergo a system restart, we'd need
      // to re-run this job.
      // On the other hand, if the DB is really down, this would be the least of our problems!
      log.error("Failed to write Job result to JobStore: {}", t)
      context stop self
      stay()
  }

  whenUnhandled {
    case Event(msg, _) =>
      log.error(s"Bad message to EngineJobExecutionActor in state $stateName($stateData): $msg")
      stay
  }

  def prepareJob(jobKey: BackendJobDescriptorKey) = {
    val jobPreparationActorName = s"$workflowId-BackendPreparationActor-${jobKey.tag}"
    val jobPrepProps = JobPreparationActor.props(executionData, jobKey, factory, initializationData, serviceRegistryActor)
    val jobPreparationActor = context.actorOf(jobPrepProps, jobPreparationActorName)
    jobPreparationActor ! JobPreparationActor.Start
    goto(PreparingJob)
  }

  private def buildJobExecutionActorName(jobDescriptor: BackendJobDescriptor) = {
    s"$workflowId-BackendJobExecutionActor-${jobDescriptor.key.tag}"
  }

  private def saveJobCompletionResult(response: BackendJobExecutionResponse) = response match {
    case SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs) => saveSuccessfulJobResults(jobKey, returnCode, jobOutputs)
    case AbortedResponse(jobKey: BackendJobDescriptorKey) => log.debug("Won't save 'aborted' job response to JobStore")
    case FailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = false)
    case FailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) => saveUnsuccessfulJobResults(jobKey, returnCode, throwable, retryable = true)
  }

  private def saveSuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], outputs: JobOutputs) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowId)
    val jobStoreResult = JobResultSuccess(returnCode, outputs)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }

  private def saveUnsuccessfulJobResults(jobKey: JobKey, returnCode: Option[Int], reason: Throwable, retryable: Boolean) = {
    val jobStoreKey = jobKey.toJobStoreKey(workflowId)
    val jobStoreResult = JobResultFailure(returnCode, reason, retryable)
    jobStoreActor ! RegisterJobCompleted(jobStoreKey, jobStoreResult)
  }
}
