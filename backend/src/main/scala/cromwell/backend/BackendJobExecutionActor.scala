package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.OutputEvaluator.EvaluatedJobOutputs
import cromwell.core.path.Path
import cromwell.core._
import wom.expression.IoFunctionSet
import wom.values.WomValue

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.util.{Success, Try}

object BackendJobExecutionActor {

  // Commands
  sealed trait BackendJobExecutionActorCommand extends BackendWorkflowLifecycleActorCommand
  case object ExecuteJobCommand extends BackendJobExecutionActorCommand
  case object RecoverJobCommand extends BackendJobExecutionActorCommand
  case object ReconnectToAbortingJobCommand extends BackendJobExecutionActorCommand
  case object ReconnectJobCommand extends BackendJobExecutionActorCommand

  // Responses
  sealed trait BackendJobExecutionActorResponse extends BackendWorkflowLifecycleActorResponse

  sealed trait BackendJobExecutionResponse extends BackendJobExecutionActorResponse { def jobKey: JobKey }
  case class JobSucceededResponse(jobKey: BackendJobDescriptorKey,
                                  returnCode: Option[Int],
                                  jobOutputs: CallOutputs,
                                  jobDetritusFiles: Option[Map[String, Path]],
                                  executionEvents: Seq[ExecutionEvent],
                                  dockerImageUsed: Option[String],
                                  resultGenerationMode: ResultGenerationMode) extends BackendJobExecutionResponse

  sealed trait ResultGenerationMode
  case object RunOnBackend extends ResultGenerationMode
  case object CallCached extends ResultGenerationMode
  case object FetchedFromJobStore extends ResultGenerationMode

  case class JobAbortedResponse(jobKey: BackendJobDescriptorKey) extends BackendJobExecutionResponse
  
  sealed trait BackendJobFailedResponse extends BackendJobExecutionResponse {  def throwable: Throwable; def returnCode: Option[Int] }
  case class JobFailedNonRetryableResponse(jobKey: JobKey, throwable: Throwable, returnCode: Option[Int]) extends BackendJobFailedResponse
  case class JobFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int], retryWithDoubleMemory: Boolean = false) extends BackendJobFailedResponse
  
  // Reconnection Exceptions
  case class JobReconnectionNotSupportedException(jobKey: BackendJobDescriptorKey) extends Exception(
    s"This backend does not support job reconnection. The status of the underlying job for ${jobKey.tag} cannot be known."
  ) with CromwellFatalExceptionMarker

  case class JobNotFoundException(jobKey: BackendJobDescriptorKey) extends Exception (
    s"No backend job for ${jobKey.tag} could be found. The status of the underlying job cannot be known."
  ) with CromwellFatalExceptionMarker

  def buildJobExecutionActorName(workflowId: WorkflowId, jobKey: BackendJobDescriptorKey) = {
    s"$workflowId-BackendJobExecutionActor-${jobKey.tag}"
  }
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendJobExecutionActor extends BackendJobLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case ExecuteJobCommand => performActionThenRespond(execute, onFailure = executionFailed)
    case RecoverJobCommand => performActionThenRespond(recover, onFailure = executionFailed)
    case ReconnectJobCommand => performActionThenRespond(reconnect, onFailure = executionFailed)
    case ReconnectToAbortingJobCommand => performActionThenRespond(reconnectToAborting, onFailure = executionFailed)
    case AbortJobCommand =>
      abort()
      context.parent ! JobAbortedResponse(jobDescriptor.key)
      context.stop(self)
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed = (t: Throwable) =>
    JobFailedNonRetryableResponse(jobKey = jobDescriptor.key, throwable = t, returnCode = None)

  /**
    * Execute a new job.
    */
  def execute: Future[BackendJobExecutionResponse]

  /**
    * Restart or resume a previously-started job.
    */
  def recover: Future[BackendJobExecutionResponse] = {
    log.warning("{} backend currently doesn't support recovering jobs. Starting {} again.",
      jobTag, jobDescriptor.key.call.fullyQualifiedName)
    execute
  }

  /**
    * Tries to reconnect to a previously started job. This method differs from recover by sending a ReconnectionFailure
    * if it can't reconnect to the job for whatever reason. It should NOT execute the job if reconnection is impossible.
    */
  def reconnect: Future[BackendJobExecutionResponse] = {
    Future.failed(JobReconnectionNotSupportedException(jobDescriptor.key))
  }

  /**
    * Similar to reconnect, except that if the reconnection succeeds and the job is still running,
    * an abort attempt should be made.
    */
  def reconnectToAborting: Future[BackendJobExecutionResponse] = {
    Future.failed(JobReconnectionNotSupportedException(jobDescriptor.key))
  }

  /**
    * Abort a running job.
    */
  def abort(): Unit = {
    log.warning("{} backend currently doesn't support abort for {}.",
      jobTag, jobDescriptor.key.call.fullyQualifiedName)
  }

  def evaluateOutputs(wdlFunctions: IoFunctionSet,
                      postMapper: WomValue => Try[WomValue] = v => Success(v))(implicit ec: ExecutionContext): EvaluatedJobOutputs = {
    Await.result(OutputEvaluator.evaluateOutputs(jobDescriptor, wdlFunctions, postMapper), Duration.Inf)
  }
}
