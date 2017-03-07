package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.actor.RobustClientHelper

import scala.concurrent.duration.FiniteDuration

trait IoClientHelper extends RobustClientHelper { this: Actor with ActorLogging with IoCommandBuilder =>
  def ioActor: ActorRef

  private [core] def ioResponseReceive: Receive = {
    case ack: IoAck[_] if hasTimeout(ack.command) =>
      cancelTimeout(ack.command)
      receive.apply(ack)
    case (context: Any, ack: IoAck[_]) if hasTimeout(context -> ack.command) =>
      cancelTimeout(context -> ack.command)
      receive.apply(context -> ack)
  }
  
  def ioReceive = robustReceive orElse ioResponseReceive
  
  def sendIoCommand(ioCommand: IoCommand[_], timeout: FiniteDuration = RobustClientHelper.DefaultRequestLostTimeout) = {
    robustSend(ioCommand, ioActor, timeout)
  }

  def sendIoCommandWithContext[T](ioCommand: IoCommand[_], context: T, timeout: FiniteDuration = RobustClientHelper.DefaultRequestLostTimeout) = {
    robustSend(context -> ioCommand, ioActor, timeout)
  }
}
