package cromwell.core.io

import akka.actor.ActorRef

/**
  * Helper trait to used async Io method easily
  */
trait AsyncIoActorClient {
  def ioActor: ActorRef
  def ioCommandBuilder: IoCommandBuilder

  val asyncIo = new AsyncIo(ioActor, ioCommandBuilder)
}
