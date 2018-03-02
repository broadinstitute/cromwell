package cromwell.services.instrumentation

import cats.data.NonEmptyVector
import cromwell.core.actor.BatchActor
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath
import cromwell.services.instrumentation.InstrumentedBatchActor.{QueueSizeTimerAction, QueueSizeTimerKey}

import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration

object InstrumentedBatchActor {
  case object QueueSizeTimerKey  
  case object QueueSizeTimerAction  
}

/**
  * Layer over batch actor that instruments the throughput and queue size
  */
abstract class InstrumentedBatchActor[C](flushRate: FiniteDuration,
                                         batchSize: Int,
                                         instrumentationPath: InstrumentationPath,
                                         instrumentationPrefix: Option[String]) extends BatchActor[C](flushRate, batchSize) 
  with CromwellInstrumentationActor {
  private val writePath = instrumentationPath.::("write")
  private val queueSizePath = instrumentationPath.::("queue")

  timers.startPeriodicTimer(QueueSizeTimerKey, QueueSizeTimerAction, CromwellInstrumentation.InstrumentationRate)
  
  whenUnhandled {
    case Event(QueueSizeTimerAction, data) =>
      sendGauge(queueSizePath, data.weight.toLong, instrumentationPrefix)
      stay()
  }
  
  protected def processInner(data: Vector[C]): Future[Int]
  
  override protected final def process(data: NonEmptyVector[C]) = {
    val action = processInner(data.toVector)
    action foreach { n =>
      count(writePath, n.toLong, instrumentationPrefix)
    }
    action
  }
}
