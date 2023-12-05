package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.io.AsyncIo.defaultTimeout
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise

import scala.concurrent.Promise
import scala.concurrent.duration.FiniteDuration

object IoPromiseProxyActor {
  case class IoCommandWithPromise[A](ioCommand: IoCommand[A], timeout: FiniteDuration = defaultTimeout) {
    val promise = Promise[A]()

    override def hashCode(): Int = ioCommand.hashCode()
  }
  def props(ioActor: ActorRef) = Props(new IoPromiseProxyActor(ioActor))
}

/**
  * Acts as a proxy for the IoActor by receiving promises along with the command and completing them when the response comes back.
  * This enables using the IoActor through promises easily from anywhere.
  * However backpressure is less efficient because the messages come back to this actor and the backpressure information can't really
  * be communicated back to the original sender (who might not be an actor)
  */
class IoPromiseProxyActor(override val ioActor: ActorRef) extends Actor with ActorLogging with IoClientHelper {
  override def receive = ioReceive orElse actorReceive

  def actorReceive: Receive = { case withPromise: IoCommandWithPromise[_] =>
    sendIoCommandWithContext(withPromise.ioCommand, withPromise.promise, withPromise.timeout)
  }

  override protected def ioResponseReceive: Receive = { case (promise: Promise[_], ack: IoAck[Any] @unchecked) =>
    cancelTimeout(promise -> ack.command)
    // This is not typesafe and assumes the Promise context is of the same type as the IoAck response.
    promise.asInstanceOf[Promise[Any]].complete(ack.toTry)
    ()
  }

  override def onTimeout(message: Any, to: ActorRef): Unit = message match {
    case (promise: Promise[_], ioAck: IoAck[_]) =>
      promise.tryFailure(IoTimeout(ioAck.command))
      ()
    case _ =>
  }
}
