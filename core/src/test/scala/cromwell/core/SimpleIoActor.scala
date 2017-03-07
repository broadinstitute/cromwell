package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.io._

import scala.io.Codec
import scala.util.{Failure, Success, Try}

object SimpleIoActor {
  def props = Props(new SimpleIoActor)
}

class SimpleIoActor extends Actor {
  
  override def receive = {
    case command: IoCopyCommand =>
      
      Try(command.source.copyTo(command.destination, command.overwrite)) match {
        case Success(_) => sender() ! IoSuccess(command, ())
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoWriteCommand =>
      
      Try(command.file.write(command.content)(command.openOptions, Codec.UTF8)) match {
        case Success(_) => sender() ! IoSuccess(command, ())
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoDeleteCommand =>
      
      Try(command.file.delete(command.swallowIOExceptions)) match {
        case Success(_) => sender() ! IoSuccess(command, ())
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoSizeCommand =>
      
      Try(command.file.size) match {
        case Success(size) => sender() ! IoSuccess(command, size)
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoContentAsStringCommand =>
      
      Try(command.file.contentAsString) match {
        case Success(content) => sender() ! IoSuccess(command, content)
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoHashCommand =>
      Try(command.file.md5) match {
        case Success(hash) => sender() ! IoSuccess(command, hash)
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }

    // With context
    case (requestContext: Any, command: IoCopyCommand) =>
      
      Try(command.source.copyTo(command.destination, command.overwrite)) match {
        case Success(_) => sender() ! (requestContext -> IoSuccess(command, ()))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
      
    case (requestContext: Any, command: IoWriteCommand) =>

      Try(command.file.write(command.content)) match {
        case Success(_) => sender() ! (requestContext -> IoSuccess(command, ()))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
      
    case (requestContext: Any, command: IoDeleteCommand) =>

      Try(command.file.delete(command.swallowIOExceptions)) match {
        case Success(_) => sender() ! (requestContext -> IoSuccess(command, ()))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
      
    case (requestContext: Any, command: IoSizeCommand) =>
      
      Try(command.file.size) match {
        case Success(size) => sender() ! (requestContext -> IoSuccess(command, size))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
      
    case (requestContext: Any, command: IoContentAsStringCommand) =>

      Try(command.file.contentAsString) match {
        case Success(content) => sender() ! (requestContext -> IoSuccess(command, content))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
      
    case (requestContext: Any, command: IoHashCommand) =>
      
      Try(command.file.md5) match {
        case Success(hash) => sender() ! (requestContext -> IoSuccess(command, hash))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }
  }
}
