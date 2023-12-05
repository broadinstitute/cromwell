package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.{Config, ConfigFactory}
import common.util.Backoff
import cromwell.core.actor.RobustClientHelper
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration._
import net.ceedubs.ficus.Ficus._

trait IoClientHelper extends RobustClientHelper { this: Actor with ActorLogging =>
  def ioActor: ActorRef

  lazy val defaultIoTimeout = RobustClientHelper.DefaultRequestLostTimeout

  protected def config = ConfigFactory.load().as[Config]("system.io.backpressure-backoff")

  override protected def initialBackoff(): Backoff = SimpleExponentialBackoff(config)

  protected def ioResponseReceive: Receive = {
    case ack: IoAck[_] if hasTimeout(ack.command) =>
      cancelTimeout(ack.command)
      receive.apply(ack)
    case (context: Any, ack: IoAck[_]) if hasTimeout(context -> ack.command) =>
      cancelTimeout(context -> ack.command)
      receive.apply(context -> ack)
  }

  def ioReceive = robustReceive orElse ioResponseReceive

  def sendIoCommand(ioCommand: IoCommand[_]) =
    sendIoCommandWithCustomTimeout(ioCommand, defaultIoTimeout)

  def sendIoCommandWithCustomTimeout(ioCommand: IoCommand[_], timeout: FiniteDuration) =
    robustSend(ioCommand, ioActor, timeout)

  def sendIoCommandWithContext[T](ioCommand: IoCommand[_], context: T, timeout: FiniteDuration = defaultIoTimeout) =
    robustSend(context -> ioCommand, ioActor, timeout)

}
