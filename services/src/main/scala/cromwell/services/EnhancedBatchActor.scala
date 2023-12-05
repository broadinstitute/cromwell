package cromwell.services

import cromwell.core.actor.BatchActor
import cromwell.services.instrumentation.{CromwellInstrumentationActor, InstrumentedBatchActor}
import cromwell.services.loadcontroller.LoadControlledBatchActor

import scala.concurrent.duration.FiniteDuration

/**
  * A BatchActor with instrumentation and load control traits mixed in to remove some boilerplate
  */
abstract class EnhancedBatchActor[C](flushRate: FiniteDuration, batchSize: Int)
    extends BatchActor[C](flushRate, batchSize)
    with InstrumentedBatchActor[C]
    with CromwellInstrumentationActor
    with LoadControlledBatchActor[C] {
  protected def enhancedReceive: Receive = loadControlReceive.orElse(instrumentationReceive).orElse(super.receive)
}
