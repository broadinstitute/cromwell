package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._
import cromwell.logging.WorkflowLogger

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success}

object CallExecutionActor {
  sealed trait CallExecutionActorMessage
  final case class IssuePollRequest(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class PollResponseReceived(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class Finish(executionHandle: ExecutionHandle) extends CallExecutionActorMessage

  sealed trait ExecutionMode extends CallExecutionActorMessage {
    def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle]
  }

  case object Execute extends ExecutionMode {
    override def execute(backendCall: BackendCall)(implicit ec: ExecutionContext) = backendCall.execute
  }

  final case class Resume(jobKey: JobKey) extends ExecutionMode {
    override def execute(backendCall: BackendCall)(implicit ec: ExecutionContext) = backendCall.resume(jobKey)
  }

  final case class UseCachedCall(cachedBackendCall: BackendCall, backendCall: BackendCall) extends ExecutionMode {
    override def execute(backendCall: BackendCall)(implicit ec: ExecutionContext) = backendCall.useCachedCall(cachedBackendCall)
  }

  def props(backendCall: BackendCall): Props = Props(new CallExecutionActor(backendCall))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(backendCall: BackendCall) extends Actor with CromwellActor {
  import CallExecutionActor._

  val akkaLogger = Logging(context.system, classOf[CallExecutionActor])
  val logger = WorkflowLogger(
    "CallExecutionActor",
    backendCall.workflowDescriptor,
    akkaLogger = Option(akkaLogger),
    callTag = Option(backendCall.key.tag)
  )

  implicit val ec = context.system.dispatcher

  /**
   * Schedule work according to the schedule of the `backoff`.
   */
  private def scheduleWork(work: => Unit): Unit = {
    val interval = backoff.nextBackOffMillis().millis
    context.system.scheduler.scheduleOnce(interval) {
      work
    }
  }

  private val backoff = new ExponentialBackOff.Builder()
    .setInitialIntervalMillis(5.seconds.toMillis.toInt)
    .setMaxIntervalMillis(30.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
    .setMultiplier(1.1)
    .build()

  /**
   * If the `work` `Future` completes successfully, perform the `onSuccess` work, otherwise schedule
   * the execution of the `onFailure` work using an exponential backoff.
   */
  def withRetry(work: Future[ExecutionHandle], onSuccess: ExecutionHandle => Unit, onFailure: => Unit): Unit = {
    work onComplete {
      case Success(s) => onSuccess(s)
      case Failure(e: GoogleJsonResponseException) if e.getStatusCode == 403 =>
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

  override def receive = LoggingReceive {
    case mode: ExecutionMode =>
      withRetry(mode.execute(backendCall),
        onSuccess = self ! IssuePollRequest(_),
        onFailure = self ! mode
      )
    case IssuePollRequest(handle) =>
      withRetry(backendCall.poll(handle),
        onSuccess = self ! PollResponseReceived(_),
        onFailure = self ! IssuePollRequest(handle)
      )
    case PollResponseReceived(handle) if handle.isDone => self ! Finish(handle)
    case PollResponseReceived(handle) => scheduleWork(self ! IssuePollRequest(handle))
    case Finish(handle) =>
      context.parent ! CallActor.ExecutionFinished(backendCall.call, handle.result)
      context.stop(self)
    case badMessage => logger.error(s"Unexpected message $badMessage.")
  }
}
