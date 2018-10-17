package cromwell.core.actor

import akka.actor.{Actor, ActorLogging, ActorRef, Cancellable}
import common.util.Backoff
import cromwell.core.actor.RobustClientHelper._
import cromwell.core.actor.StreamIntegration._
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration.{FiniteDuration, _}
import scala.language.postfixOps

object RobustClientHelper {
  case class RequestTimeout(msg: Any, to: ActorRef)
  val DefaultRequestLostTimeout = 5 minutes
}

trait RobustClientHelper { this: Actor with ActorLogging =>
  private [actor] implicit val robustActorHelperEc = context.dispatcher

  private var backoff: Option[Backoff] = None

  // package private for testing
  private [core] var timeouts = Map.empty[Any, (Cancellable, FiniteDuration)]

  protected def initialBackoff(): Backoff = SimpleExponentialBackoff(5.seconds, 20.minutes, 2D)

  def robustReceive: Receive = {
    case BackPressure(request) =>
      val snd = sender()
      newTimer(request, snd, generateBackpressureTime)
      resetTimeout(request, snd)
      ()
    case RequestTimeout(request, to) => onTimeout(request, to)
  }

  private final def newTimer(msg: Any, to: ActorRef, in: FiniteDuration) = {
    context.system.scheduler.scheduleOnce(in, to, msg)(robustActorHelperEc, self)
  }

  def robustSend(msg: Any, to: ActorRef, timeout: FiniteDuration = DefaultRequestLostTimeout): Unit = {
    to ! msg
    addTimeout(msg, to, timeout)
  }

  private final def addTimeout(command: Any, to: ActorRef, timeout: FiniteDuration) = {
    val cancellable = newTimer(RequestTimeout(command, to), self, timeout)
    timeouts = timeouts + (command -> (cancellable -> timeout))
  }

  protected final def hasTimeout(command: Any) = timeouts.get(command).isDefined

  protected final def cancelTimeout(command: Any) = {
    timeouts.get(command) foreach { case (cancellable, _) => cancellable.cancel() }
    timeouts = timeouts - command
  }

  private final def resetTimeout(command: Any, to: ActorRef) = {
    val timeout = timeouts.get(command) map { _._2 }
    cancelTimeout(command)
    timeout foreach { addTimeout(command, to, _) }
  }

  private [actor] final def generateBackpressureTime: FiniteDuration = {
    val effectiveBackoff = backoff.getOrElse({
      val firstBackoff = initialBackoff()
      backoff = Option(firstBackoff)
      firstBackoff
    })
    val backoffTime = effectiveBackoff.backoffMillis
    backoff = Option(effectiveBackoff.next)
    backoffTime.millis
  }

  protected def onTimeout(message: Any, to: ActorRef): Unit
}
