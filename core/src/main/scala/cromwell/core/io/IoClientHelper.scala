package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.actor.RobustClientHelper

import scala.concurrent.duration.FiniteDuration

trait IoClientHelper extends RobustClientHelper { this: Actor with ActorLogging with IoCommandBuilder =>
  def ioActor: ActorRef
  
  def sendIoCommand(ioCommand: IoCommand[_], timeout: FiniteDuration = RobustClientHelper.DefaultRequestLostTimeout) = {
    robustSend(ioCommand, ioActor, timeout)
  }

  def sendIoCommandWithContext[T](ioCommand: IoCommand[_], context: T, timeout: FiniteDuration = RobustClientHelper.DefaultRequestLostTimeout) = {
    robustSend(context -> ioCommand, ioActor, timeout)
  }
}
