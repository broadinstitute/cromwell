package cromwell.server

import akka.event.DeadLetterListener

class CromwellDeadLetterListener extends DeadLetterListener {
  def shutdownReceive: Receive = {
    // This silences the dead letter messages when Cromwell is shutting down
    case _ if CromwellShutdown.shutdownInProgress() =>
  }
  override def receive = shutdownReceive.orElse(super.receive)
}
