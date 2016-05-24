package cromwell.engine.callexecution

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._
import cromwell.engine.callactor.CallActor
import cromwell.engine.callexecution.CallExecutionActor.{ExecutionMode, Finish, IssuePollRequest, PollResponseReceived}
import cromwell.engine.finalcall.FinalCall
import cromwell.engine.{CromwellActor, CromwellFatalException}
import cromwell.logging.WorkflowLogger
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s.Scope

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object CallExecutionActor {
  sealed trait CallExecutionActorMessage
  final case class IssuePollRequest(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class PollResponseReceived(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class Finish(executionHandle: ExecutionHandle) extends CallExecutionActorMessage

  sealed trait ExecutionMode extends CallExecutionActorMessage

  case object Execute extends ExecutionMode
  final case class Resume(jobKey: BackendJobKey) extends ExecutionMode
  final case class UseCachedCall(cachedBackendCall: BackendCallJobDescriptor) extends ExecutionMode

  def props(backendCall: BackendCallJobDescriptor): Props = Props(new BackendCallExecutionActor(backendCall))
  def props(finalCall: FinalCall, workflowMetadataResponse: WorkflowMetadataResponse): Props = {
    Props(new FinalCallExecutionActor(finalCall, workflowMetadataResponse))
  }
}

trait CallExecutionActor extends Actor with CromwellActor {

  val akkaLogger = Logging(context.system, classOf[CallExecutionActor])
  def logger: WorkflowLogger

  implicit val ec = context.system.dispatcher

  /**
    * Schedule work according to the schedule of the `backoff`.
    */
  protected def scheduleWork(work: => Unit): Unit = {
    val interval = backoff.nextBackOffMillis().millis
    context.system.scheduler.scheduleOnce(interval) {
      work
    }
  }

  def backoff: ExponentialBackOff

  /**
    * If the `work` `Future` completes successfully, perform the `onSuccess` work, otherwise schedule
    * the execution of the `onFailure` work using an exponential backoff.
    */
  def withRetry(work: Future[ExecutionHandle], onSuccess: ExecutionHandle => Unit, onFailure: => Unit): Unit = {
    work onComplete {
      case Success(s) => onSuccess(s)
      case Failure(e: CromwellFatalException) =>
        logger.error(e.getMessage, e)
        self ! Finish(FailedExecutionHandle(e))
      case Failure(e: Exception) =>
        logger.error(e.getMessage, e)
        scheduleWork(onFailure)
      case Failure(throwable) =>
        // This is a catch-all for a JVM-ending kind of exception, which is why we throw the exception
        logger.error(throwable.getMessage, throwable)
        throw throwable
    }
  }

  /**
    * Update the ExecutionHandle
    */
  def poll(handle: ExecutionHandle): Future[ExecutionHandle]

  /**
    * Start the execution. Once the Future resolves, the ExecutionHandle can be used to poll
    * the state of the execution.
    */
  def execute(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle]
  def call: Scope

  override def receive = LoggingReceive {
    case mode: ExecutionMode =>
      withRetry(execute(mode),
        onSuccess = self ! IssuePollRequest(_),
        onFailure = self ! mode
      )

    case IssuePollRequest(handle) =>
      withRetry(poll(handle),
        onSuccess = self ! PollResponseReceived(_),
        onFailure = self ! IssuePollRequest(handle)
      )
    case PollResponseReceived(handle) if handle.isDone => self ! Finish(handle)
    case PollResponseReceived(handle) => scheduleWork(self ! IssuePollRequest(handle))
    case Finish(handle) =>
      context.parent ! CallActor.ExecutionFinished(call, handle.result)
      context.stop(self)
    case badMessage => logger.error(s"Unexpected message $badMessage.")
  }
}