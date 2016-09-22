package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.{ExecutionEvent, JobOutputs}
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

  sealed trait BackendJobExecutionResponse extends BackendJobExecutionActorResponse { def jobKey: BackendJobDescriptorKey }
  case class SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs, jobDetritusFiles: Option[Map[String, String]], executionEvents: Seq[ExecutionEvent]) extends BackendJobExecutionResponse
  case class AbortedResponse(jobKey: BackendJobDescriptorKey) extends BackendJobExecutionResponse
  case class FailedNonRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) extends BackendJobExecutionResponse
  case class FailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable, returnCode: Option[Int]) extends BackendJobExecutionResponse
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
    FailedNonRetryableResponse(jobKey = jobDescriptor.key, throwable = t, returnCode = None)

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
