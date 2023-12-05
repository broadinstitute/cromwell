package cromwell.engine.io.gcs

import akka.actor.Scheduler
import akka.stream.FlowShape
import akka.stream.scaladsl.{Balance, GraphDSL, Merge}
import cromwell.engine.io.IoActor.IoResult
import cromwell.engine.io.IoCommandContext
import cromwell.engine.io.gcs.GcsBatchFlow.GcsBatchFlowConfig

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.FiniteDuration

/**
  * Balancer that distributes requests to multiple batch flows in parallel
  */
class ParallelGcsBatchFlow(config: GcsBatchFlowConfig,
                           scheduler: Scheduler,
                           onRetry: IoCommandContext[_] => Throwable => Unit,
                           onBackpressure: Option[Double] => Unit,
                           applicationName: String,
                           commandBackpressureStaleness: FiniteDuration
)(implicit ec: ExecutionContext) {

  // noinspection TypeAnnotation
  val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._
    val balancer = builder.add(Balance[GcsBatchCommandContext[_, _]](config.parallelism, waitForAllDownstreams = false))
    val merge = builder.add(Merge[IoResult](config.parallelism))

    for (_ <- 1 to config.parallelism) {
      val workerFlow = new GcsBatchFlow(
        batchSize = config.maxBatchSize,
        batchTimespan = config.maxBatchDuration,
        scheduler = scheduler,
        onRetry = onRetry,
        onBackpressure = onBackpressure,
        applicationName = applicationName,
        backpressureStaleness = commandBackpressureStaleness
      ).flow
      // for each worker, add an edge from the balancer to the worker, then wire
      // it to the merge element
      balancer ~> workerFlow.async ~> merge
    }

    FlowShape(balancer.in, merge.out)
  }
}
