package cromwell.backend

import akka.actor.ActorLogging
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, _}
import cromwell.backend.BackendLifecycleActor._
import cromwell.core.{CallOutput, CallOutputs}
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
  case class BackendJobExecutionSucceededResponse(jobKey: BackendJobDescriptorKey, callOutputs: CallOutputs) extends BackendJobExecutionResponse
  case class BackendJobExecutionAbortedResponse(jobKey: BackendJobDescriptorKey) extends BackendJobExecutionResponse
  case class BackendJobExecutionFailedResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends BackendJobExecutionResponse
  case class BackendJobExecutionFailedRetryableResponse(jobKey: BackendJobDescriptorKey, throwable: Throwable) extends BackendJobExecutionResponse
}

/**
  * Workflow-level actor for executing, recovering and aborting jobs.
  */
trait BackendJobExecutionActor extends BackendJobLifecycleActor with ActorLogging {

  def receive: Receive = LoggingReceive {
    case ExecuteJobCommand => performActionThenRespond(execute, onFailure = executionFailed)
    case RecoverJobCommand => performActionThenRespond(recover, onFailure = executionFailed)
    case AbortJobCommand => abortJob
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed = (t: Throwable) => BackendJobExecutionFailedResponse(jobDescriptor.key, t)

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
  def abortJob: Unit

  def evaluateOutputs(wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)) = {
    val inputs = jobDescriptor.inputs
    jobDescriptor.call.task.outputs.foldLeft(Map.empty[LocallyQualifiedName, Try[CallOutput]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (name, value) if value.isSuccess => name -> value.get.wdlValue
      }
      def lookup = (currentOutputs ++ inputs).apply _
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val callOutput = output.name -> (coerced flatMap postMapper map { CallOutput(_, None) })

      outputMap + callOutput

    })
  }
}
