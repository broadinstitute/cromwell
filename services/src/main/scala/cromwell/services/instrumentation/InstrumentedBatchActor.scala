package cromwell.services.instrumentation

import cats.data.NonEmptyList
import cromwell.core.actor.BatchActor
import cromwell.services.instrumentation.InstrumentedBatchActor.{QueueSizeTimerAction, QueueSizeTimerKey}

import scala.concurrent.Future
import scala.util.{Failure, Success}

object InstrumentedBatchActor {
  case object QueueSizeTimerKey
  case object QueueSizeTimerAction
}

/**
  * Can be mixed in a BatchActor to provide helper methods to instrument queue size and throughput
  */
trait InstrumentedBatchActor[C] { this: BatchActor[C] with CromwellInstrumentation =>

  protected def instrumentationPath: NonEmptyList[String]
  protected def instrumentationPrefix: Option[String]

  // If this actor is behind a router, add its name to the instrumentation path so that all routees don't override each other's values
  private def makePath(name: String) = if (routed)
    instrumentationPath.concatNel(NonEmptyList.of(self.path.name, name))
  else
    instrumentationPath.concatNel(NonEmptyList.one(name))

  private val processedPath = makePath("processed")
  private val failurePath = makePath("failure")
  private val queueSizePath = makePath("queue")

  timers.startSingleTimer(QueueSizeTimerKey, QueueSizeTimerAction, CromwellInstrumentation.InstrumentationRate)

  /**
    * Don't forget to chain this into your receive method to instrument the queue size:
    * override def receive = instrumentationReceive.orElse(super.receive)
    */
  protected def instrumentationReceive: Receive = {
    case QueueSizeTimerAction => 
      sendGauge(queueSizePath, stateData.weight.toLong, instrumentationPrefix)
      timers.startSingleTimer(QueueSizeTimerKey, QueueSizeTimerAction, CromwellInstrumentation.InstrumentationRate)
  }

  /**
    * Don't forget to wrap your `process` or `processHead` method with this function if you want
    * to instrument your processing rate:
    * instrumentedProcess {
    *   do work
    * }
    */
  protected def instrumentedProcess(f: => Future[Int]) = {
    val action = f
    action onComplete {
      case Success(numProcessed) => count(processedPath, numProcessed.toLong, instrumentationPrefix)
      case Failure(_) => count(failurePath, 1L, instrumentationPrefix)
    }
    action
  }
}
