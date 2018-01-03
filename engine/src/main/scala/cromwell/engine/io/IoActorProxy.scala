package cromwell.engine.io

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.io.{IoCommand, IoPromiseProxyActor}
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

object IoActorProxy {
  def props(ioActor: ActorRef) = Props(new IoActorProxy(ioActor))
}

class IoActorProxy(ioActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper {
  private val ioPromiseProxyActor: ActorRef = context.actorOf(IoPromiseProxyActor.props(ioActor), "IoPromiseProxyActor")

  override def receive = {
    // If it's an IoCommandWithPromise send it to the proxy actor
    case ioPromise: IoCommandWithPromise[_] => ioPromiseProxyActor forward ioPromise
    // otherwise to the IoActor directly
    case ioCommand: IoCommand[_] => ioActor forward ioCommand
    case withContext: (Any, IoCommand[_]) @unchecked => ioActor forward withContext

    case ShutdownCommand => 
      context stop ioPromiseProxyActor
      waitForActorsAndShutdown(NonEmptyList.one(ioActor))
  }
}
