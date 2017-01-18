package cromwell.services.io

import akka.actor.Actor

abstract class IoServiceActor extends Actor {

  override def receive = {
    case command: IoActorCommand[_] => handleCommand(command) foreach { sender ! _ }
  }

  def handleCommand(command: IoActorCommand[_]): Option[Any]
}
