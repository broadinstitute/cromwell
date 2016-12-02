package cromwell.backend

import java.nio.file.Path

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.wdl.OutputEvaluator
import cromwell.core.{CallOutputs, ExecutionEvent, JobKey}
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.util.{Success, Try}

object BackendJobExecutionActor {

  // Commands
  sealed trait BackendJobExecutionActorCommand extends BackendWorkflowLifecycleActorCommand
  case object ExecuteJobCommand extends BackendJobExecutionActorCommand
  case object RecoverJobCommand extends BackendJobExecutionActorCommand

  // Responses
  sealed trait BackendJobExecutionActorResponse extends BackendWorkflowLifecycleActorResponse

  sealed trait BackendJobExecutionResponse extends BackendJobExecutionActorResponse { def jobKey: JobKey }
  case class JobSucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: CallOutputs, jobDetritusFiles: Option[Map[String, Path]], executionEvents: Seq[ExecutionEvent]) extends BackendJobExecutionResponse
  case class AbortedResponse(jobKey: BackendJobDescriptorKey) extends BackendJobExecutionResponse
  sealed trait BackendJobFailedResponse extends BackendJobExecutionResponse {  def throwable: Throwable; def returnCode: Option[Int] }
  case class JobFailedNonRetryableResponse(jobKey: JobKey, throwable: Throwable, returnCode: Option[Int]) extends BackendJobFailedResponse
  case class JobFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) extends BackendJobFailedResponse
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendJobExecutionActor extends BackendJobLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case ExecuteJobCommand => performActionThenRespond(execute, onFailure = executionFailed)
    case RecoverJobCommand => performActionThenRespond(recover, onFailure = executionFailed)
    case AbortJobCommand =>
      abort()
      context.parent ! AbortedResponse(jobDescriptor.key)
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
    * Abort a running job.
    */
  def abort(): Unit = {
    log.warning("{} backend currently doesn't support abort for {}.",
      jobTag, jobDescriptor.key.call.fullyQualifiedName)
  }

  def evaluateOutputs(wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)) = {
    OutputEvaluator.evaluateOutputs(jobDescriptor, wdlFunctions, postMapper)
  }
}
