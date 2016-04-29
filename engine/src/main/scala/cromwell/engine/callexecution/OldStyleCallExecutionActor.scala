package cromwell.engine.callexecution

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._
import cromwell.engine.callactor.OldStyleCallActor
import cromwell.engine.callexecution.OldStyleCallExecutionActor.{ExecutionMode, Finish, IssuePollRequest, PollResponseReceived}
import cromwell.engine.finalcall.OldStyleFinalCall
import cromwell.engine.{CromwellActor, CromwellFatalException}
import cromwell.logging.WorkflowLogger
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s.Scope

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleCallExecutionActor {
  sealed trait CallExecutionActorMessage
  final case class IssuePollRequest(executionHandle: OldStyleExecutionHandle) extends CallExecutionActorMessage
  final case class PollResponseReceived(executionHandle: OldStyleExecutionHandle) extends CallExecutionActorMessage
  final case class Finish(executionHandle: OldStyleExecutionHandle) extends CallExecutionActorMessage

  sealed trait ExecutionMode extends CallExecutionActorMessage

  case object Execute extends ExecutionMode
  final case class Resume(executionInfos: Map[String, Option[String]]) extends ExecutionMode
  final case class UseCachedCall(cachedBackendCall: OldStyleBackendCallJobDescriptor) extends ExecutionMode

  def props(backendCall: OldStyleBackendCallJobDescriptor): Props = Props(new OldStyleBackendCallExecutionActor(backendCall))
  def props(finalCall: OldStyleFinalCall, workflowMetadataResponse: WorkflowMetadataResponse): Props = {
    Props(new FinalCallExecutionActor(finalCall, workflowMetadataResponse))
  }
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait OldStyleCallExecutionActor extends Actor with CromwellActor {

  val akkaLogger = Logging(context.system, classOf[OldStyleCallExecutionActor])
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
  def withRetry(work: Future[OldStyleExecutionHandle], onSuccess: OldStyleExecutionHandle => Unit, onFailure: => Unit): Unit = {
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
  def poll(handle: OldStyleExecutionHandle): Future[OldStyleExecutionHandle]

  /**
    * Start the execution. Once the Future resolves, the ExecutionHandle can be used to poll
    * the state of the execution.
    */
  def execute(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle]
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
      context.parent ! OldStyleCallActor.ExecutionFinished(call, handle.result)
      context.stop(self)
    case badMessage => logger.error(s"Unexpected message $badMessage.")
  }
}