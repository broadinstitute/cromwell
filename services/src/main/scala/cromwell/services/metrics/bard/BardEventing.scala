package cromwell.services.metrics.bard

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.metrics.bard.model.BardEvent

object BardEventing {
  sealed trait BardEventingMessage extends ServiceRegistryMessage {
    override def serviceName: String = "BardEventing"
  }

  case class BardEventRequest(event: BardEvent) extends BardEventingMessage

}
