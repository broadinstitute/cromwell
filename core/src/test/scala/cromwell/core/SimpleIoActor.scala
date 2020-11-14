package cromwell.core

import java.nio.charset.StandardCharsets

import akka.actor.{Actor, Props}
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.io._

import scala.concurrent.Promise
import scala.util.{Failure, Success, Try}

object SimpleIoActor {
  def props: Props = Props(new SimpleIoActor)
}

class SimpleIoActor extends Actor {
  
  override def receive: Receive = {
    case command: IoCopyCommand =>
      
      Try(command.source.copyTo(command.destination)) match {
        case Success(_) => sender() ! IoSuccess(command, ())
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }
      
    case command: IoWriteCommand =>
      
      Try(command.file.write(command.content)(command.openOptions, StandardCharsets.UTF_8)) match {
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
      
    case command: IoExistsCommand =>
      Try(command.file.exists) match {
        case Success(exists) => sender() ! IoSuccess(command, exists)
        case Failure(failure) => sender() ! IoFailure(command, failure)
      }

    // With context
    case (requestContext: Any, command: IoCopyCommand) =>
      
      Try(command.source.copyTo(command.destination, overwrite = true)) match {
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

    case (requestContext: Any, command: IoExistsCommand) =>

      Try(command.file.exists) match {
        case Success(exists) => sender() ! (requestContext -> IoSuccess(command, exists))
        case Failure(failure) => sender() ! (requestContext -> IoFailure(command, failure))
      }

    case withPromise: IoCommandWithPromise[_] => self ! ((withPromise.promise, withPromise.ioCommand))

    case (promise: Promise[_], ack: IoAck[Any] @unchecked) =>
      promise.asInstanceOf[Promise[Any]].complete(ack.toTry)
      ()
  }
}
