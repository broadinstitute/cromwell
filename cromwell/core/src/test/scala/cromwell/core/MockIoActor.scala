package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.io._

object MockIoActor {
  def props(returnCode: String = "0", stderrSize: Long = 0L) = Props(new MockIoActor(returnCode, stderrSize))
}

class MockIoActor(returnCode: String, stderrSize: Long) extends Actor {
  override def receive = {
    case command: IoCopyCommand => sender() ! IoSuccess(command, ())
    case command: IoWriteCommand => sender() ! IoSuccess(command, ())
    case command: IoDeleteCommand => sender() ! IoSuccess(command, ())
    case command: IoSizeCommand => sender() ! IoSuccess(command, 0L)
    case command: IoContentAsStringCommand => sender() ! IoSuccess(command, "0")
      
      // With context
    case (requestContext: Any, command: IoCopyCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoWriteCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoDeleteCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoSizeCommand) => sender() ! (requestContext -> IoSuccess(command, stderrSize))
    case (requestContext: Any, command: IoContentAsStringCommand) => sender() ! (requestContext -> IoSuccess(command, returnCode))
  }
}
