package cromwell.core

import akka.actor.{Actor, Props}
import cromwell.core.FailIoActor._
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.io._

import scala.concurrent.Promise

object FailIoActor {
  def props() = Props(new FailIoActor())
  val failure = new Exception("Io failure - part of test flow")
}

class FailIoActor() extends Actor {
  override def receive = {
    case command: IoCommand[_] => sender() ! IoFailure(command, failure)

    // With context
    case (requestContext: Any, command: IoCommand[_]) => sender() ! (requestContext -> IoFailure(command, failure))
    case withPromise: IoCommandWithPromise[_] => self ! ((withPromise.promise, withPromise.ioCommand))
    case (promise: Promise[_], ack: IoAck[Any] @unchecked) =>
      promise.asInstanceOf[Promise[Any]].complete(ack.toTry)
      ()
  }
}
