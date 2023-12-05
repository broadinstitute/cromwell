package cromwell.services

import cromwell.core.actor.ThrottlerActor
import cromwell.services.instrumentation.{CromwellInstrumentationActor, InstrumentedBatchActor}
import cromwell.services.loadcontroller.LoadControlledBatchActor

/**
  * A ThrottlerActor with instrumentation and load control traits mixed in to remove some boilerplate
  */
abstract class EnhancedThrottlerActor[C]
    extends ThrottlerActor[C]
    with InstrumentedBatchActor[C]
    with CromwellInstrumentationActor
    with LoadControlledBatchActor[C] {
  protected def enhancedReceive: Receive = loadControlReceive.orElse(instrumentationReceive).orElse(super.receive)
}
