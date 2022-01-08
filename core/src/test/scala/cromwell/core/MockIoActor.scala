package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.io._

import scala.concurrent.Promise

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
    case command: IoExistsCommand => sender() ! IoSuccess(command, false)
      
      // With context
    case (requestContext: Any, command: IoCopyCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoWriteCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoDeleteCommand) => sender() ! (requestContext -> IoSuccess(command, ()))
    case (requestContext: Any, command: IoSizeCommand) => sender() ! (requestContext -> IoSuccess(command, stderrSize))
    case (requestContext: Any, command: IoContentAsStringCommand) => sender() ! (requestContext -> IoSuccess(command, returnCode))
    case (requestContext: Any, command: IoExistsCommand) => sender() ! (requestContext -> IoSuccess(command, false))

    case withPromise: IoCommandWithPromise[_] => self ! ((withPromise.promise, withPromise.ioCommand))

    case (promise: Promise[_], ack: IoAck[Any] @unchecked) =>
      promise.asInstanceOf[Promise[Any]].complete(ack.toTry)
      ()
  }
}
