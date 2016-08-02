package cromwell.backend

import akka.actor.{ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend.BackendLifecycleActor._
import cromwell.backend.callcaching.{CallCachingActor, SyncCallCachingActor}
import cromwell.backend.callcaching.CallCachingActor._
import cromwell.core.JobOutputs
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

  val callCacheActor = context.actorOf(SyncCallCachingActor.props(jobDescriptor, _ => "Implement me!", _ => "Implement me!", readFromCache = false, writeToCache = false))

  private var executeClient: Option[ActorRef] = None

  def receive: Receive = LoggingReceive {
    case ExecuteJobCommand =>
      executeClient = Option(sender)
      callCacheActor ! CheckCache
    case CallCacheHit(backendJobExecutionResponse: BackendJobExecutionResponse) => executeClient foreach { _ ! backendJobExecutionResponse }
    case CallCacheMiss =>
      log.info(s"Call Cache miss for ${jobDescriptor.call.fullyQualifiedName}:${jobDescriptor.key.index}:${jobDescriptor.key.attempt}")
      performActionThenRespond(executeWithCacheWriteOn, onFailure = executionFailed, respondTo = executeClient.get)

    case RecoverJobCommand => performActionThenRespond(recover, onFailure = executionFailed)
    case AbortJobCommand =>
      abort()
      context.parent ! AbortedResponse(jobDescriptor.key)
      context.stop(self)
  }

  // We need this for receive because we can't do `onFailure = ExecutionFailure` directly - because BackendJobDescriptor =/= BackendJobDescriptorKey
  private def executionFailed = (t: Throwable) =>
    FailedNonRetryableResponse(jobKey = jobDescriptor.key, throwable = t, returnCode = None)

  private def executeWithCacheWriteOn = execute map {
    case x @ SucceededResponse(jobKey: BackendJobDescriptorKey, returnCode: Option[Int], jobOutputs: JobOutputs) =>
      callCacheActor ! WriteCache(returnCode, jobOutputs)
      x
    case other => other
  }

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
    log.warning("{} backend currently doesn't support abort.",
      jobTag, jobDescriptor.key.call.fullyQualifiedName)
  }

  def evaluateOutputs(wdlFunctions: WdlStandardLibraryFunctions,
                      postMapper: WdlValue => Try[WdlValue] = v => Success(v)) = {
    OutputEvaluator.evaluateOutputs(jobDescriptor, wdlFunctions, postMapper)
  }
}
