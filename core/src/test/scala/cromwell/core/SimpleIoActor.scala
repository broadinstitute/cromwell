package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.io._

object SimpleIoActor {
  def props = Props(new SimpleIoActor)
}

class SimpleIoActor extends Actor {
  override def receive = {
    case command: IoCopyCommand =>
      command.source.copyTo(command.destination)
      sender() ! IoSuccess(command, ())
    case command: IoWriteCommand =>
      command.file.write(command.content)
      sender() ! IoSuccess(command, ())
    case command: IoDeleteCommand =>
      command.file.delete()
      sender() ! IoSuccess(command, ())
    case command: IoSizeCommand =>
      sender() ! IoSuccess(command, command.file.size)
    case command: IoContentAsStringCommand =>
      sender() ! IoSuccess(command, command.file.contentAsString)

    // With context
    case (requestContext: Any, command: IoCopyCommand) =>
      command.source.copyTo(command.destination)
      sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoWriteCommand) =>
      command.file.write(command.content)
      sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoDeleteCommand) =>
      command.file.delete()
      sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoSizeCommand) => 
      sender() ! (requestContext -> IoSuccess(command, command.file.size))
    case (requestContext: Any, command: IoContentAsStringCommand) => 
      sender() ! (requestContext -> IoSuccess(command, command.file.contentAsString))
  }
}
