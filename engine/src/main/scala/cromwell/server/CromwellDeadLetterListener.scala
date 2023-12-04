package cromwell.server

import akka.actor.{ActorLogging, ActorRef, DeadLetter}
import akka.event.DeadLetterListener
import cats.Show
import cats.syntax.show._

class CromwellDeadLetterListener extends DeadLetterListener with ActorLogging {
  implicit val showActor: Show[ActorRef] =
    Show.show(actor => s"Actor of path ${actor.path} toString ${actor.toString()}")

  def shutdownReceive: Receive = {
    // This silences the dead letter messages when Cromwell is shutting down
    case DeadLetter(msg, from, to) if CromwellShutdown.shutdownInProgress() =>
      log.debug(
        s"Got a dead letter during Cromwell shutdown. Sent by\n${from.show}\nto ${to.show}\n consisting of message: $msg\n "
      )
  }
  override def receive = shutdownReceive.orElse(super.receive)
}
