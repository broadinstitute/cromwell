package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.FailIoActor._
import cromwell.core.io._

object FailIoActor {
  def props() = Props(new FailIoActor())
  val failure = new Exception("Io failure - part of test flow")
}

class FailIoActor() extends Actor {
  override def receive = {
    case command: IoCommand[_] => sender() ! IoFailure(command, failure)

    // With context
    case (requestContext: Any, command: IoCommand[_]) => sender() ! (requestContext -> IoFailure(command, failure))
  }
}
