package cromwell.engine.io.gcs

import akka.actor.Scheduler
import akka.stream.FlowShape
import akka.stream.scaladsl.{Balance, GraphDSL, Merge}
import cromwell.engine.io.IoActor.IoResult
import cromwell.engine.io.IoCommandContext

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.FiniteDuration

/**
  * Balancer that distributes requests to multiple batch flows in parallel
  */
class ParallelGcsBatchFlow(parallelism: Int,
                           batchSize: Int,
                           batchTimespan: FiniteDuration,
                           scheduler: Scheduler,
                           onRetry: IoCommandContext[_] => Throwable => Unit,
                           onBackpressure: Option[Double] => Unit,
                           applicationName: String)
                          (implicit ec: ExecutionContext) {

  //noinspection TypeAnnotation
  val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._
    val balancer = builder.add(Balance[GcsBatchCommandContext[_, _]](parallelism, waitForAllDownstreams = false))
    val merge = builder.add(Merge[IoResult](parallelism))

    for (_ <- 1 to parallelism) {
      val workerFlow = new GcsBatchFlow(batchSize, batchTimespan, scheduler, onRetry, onBackpressure, applicationName).flow
      // for each worker, add an edge from the balancer to the worker, then wire
      // it to the merge element
      balancer ~> workerFlow.async ~> merge
    }

    FlowShape(balancer.in, merge.out)
  }
}
