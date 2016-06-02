package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor.{FailedNonRetryableResponse, _}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.{JobOutput, JobOutputs}
import wdl4s._
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

  sealed trait BackendJobExecutionResponse extends BackendJobExecutionActorResponse
  case class SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs) extends BackendJobExecutionResponse
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
      abort
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
  def recover: Future[BackendJobExecutionResponse]

  /**
    * Abort a running job.
    */
  def abort: Unit

  def evaluateOutputs(wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)) = {
    val inputs = jobDescriptor.inputs
    jobDescriptor.call.task.outputs.foldLeft(Map.empty[LocallyQualifiedName, Try[JobOutput]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (name, value) if value.isSuccess => name -> value.get.wdlValue
      }
      def lookup = (currentOutputs ++ inputs).apply _
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val jobOutput = output.name -> (coerced flatMap postMapper map { JobOutput(_, None) })

      outputMap + jobOutput

    })
  }
}
