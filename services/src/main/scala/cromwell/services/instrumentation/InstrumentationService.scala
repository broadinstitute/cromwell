package cromwell.services.instrumentation

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

object InstrumentationService {
  case class InstrumentationServiceMessage(metric: CromwellMetric) extends ServiceRegistryMessage {
    override val serviceName = "Instrumentation"
  }
}
