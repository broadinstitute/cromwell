package cromwell.backend.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem, Props}
import cromwell.backend.callcaching.BackendSpecificHasherActor.SingleFileHashRequest
import cromwell.core.callcaching.HashingFailedMessage

trait BackendHashingMethods {
  def hashableRuntimeAttributes: List[String]
  val fileContentsHasherActor: ActorRef
}

case class DefaultBackendHashingMethods(actorSystem: ActorSystem) extends BackendHashingMethods {
  override def hashableRuntimeAttributes = List.empty
  override val fileContentsHasherActor: ActorRef = actorSystem.actorOf(Props(new DefaultFileContentsHasherActor))
}

class DefaultFileContentsHasherActor extends Actor with ActorLogging {
  override def receive = {
    case x: SingleFileHashRequest =>
      sender ! HashingFailedMessage(x.hashKey, new NotImplementedError("This backend does not support file hashing"))
  }
}