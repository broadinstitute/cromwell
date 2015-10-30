package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object CallExecutionActor {
  sealed trait CallExecutionActorMessage
  final case class IssuePollRequest(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class PollResponseReceived(executionHandle: ExecutionHandle) extends CallExecutionActorMessage
  final case class Finish(executionHandle: ExecutionHandle) extends CallExecutionActorMessage

  sealed trait ExecutionMode extends CallExecutionActorMessage {
    def execute(backendCall: BackendCall): Future[ExecutionHandle]
  }

  case object Execute extends ExecutionMode {
    override def execute(backendCall: BackendCall) = backendCall.execute
  }

  final case class Resume(jobKey: JobKey) extends ExecutionMode {
    override def execute(backendCall: BackendCall) = backendCall.resume(jobKey)
  }

  def props(backendCall: BackendCall): Props = Props(new CallExecutionActor(backendCall))
}

/** Actor to manage the execution of a single call. */
class CallExecutionActor(backendCall: BackendCall) extends Actor with CromwellActor {
  import CallExecutionActor._

  private val log = Logging(context.system, classOf[CallExecutionActor])

  implicit val ec = context.system.dispatcher

  /**
   * Schedule an `IssuePollRequest` message parameterized by `handle` to be delivered to this actor according to
   * the schedule of the `backoff`.
   */
  private def scheduleNextPoll(handle: ExecutionHandle): Unit = {
    val interval = backoff.nextBackOffMillis().millis
    context.system.scheduler.scheduleOnce(interval) {
      self ! IssuePollRequest(handle)
    }
  }

  private val backoff = new ExponentialBackOff.Builder()
    .setInitialIntervalMillis(5.seconds.toMillis.toInt)
    .setMaxIntervalMillis(1.minute.toMillis.toInt)
    .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
    .setMultiplier(1.1)
    .build()

  val tag = s"CallExecutionActor [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.key.tag}]"

  def withHandle(successCode: ExecutionHandle => Unit): Try[ExecutionHandle] => Unit = {
    case Success(handle) => successCode(handle)
    case Failure(t) =>
      log.error(t.getMessage, t)
      self ! Finish(FailedExecutionHandle(t))
  }

  override def receive = LoggingReceive {
    case mode: ExecutionMode =>
      mode.execute(backendCall) onComplete withHandle { handle =>
        self ! IssuePollRequest(handle)
      }
    case IssuePollRequest(handle) =>
      backendCall.poll(handle) onComplete withHandle { nextHandle =>
        self ! PollResponseReceived(nextHandle)
      }
    case PollResponseReceived(handle) if handle.isDone =>
      self ! Finish(handle)
    case PollResponseReceived(handle) =>
      scheduleNextPoll(handle)
    case Finish(handle) =>
      context.parent ! CallActor.ExecutionFinished(backendCall.call, handle.result)
      context.stop(self)
    case badMessage => log.error(s"$tag: unexpected message $badMessage.")
  }
}
