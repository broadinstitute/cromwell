package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{apply => _, _}
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.util.SampleWdl
import org.scalameter.api._
import wdl4s.{TaskCall, WdlNamespaceWithWorkflow}
import org.scalameter.picklers.Implicits._

/**
  * Benchmarks the performance of the execution store.
  */
object ExecutionStoreBenchmark extends Bench[Double] {

  /* Benchmark configuration */
  lazy val measurer = new Measurer.Default
  lazy val executor = SeparateJvmsExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  lazy val reporter = new LoggingReporter[Double]
  lazy val persistor = Persistor.None
  
  val wdl = WdlNamespaceWithWorkflow.load(SampleWdl.PrepareScatterGatherWdl().wdlSource(), Seq.empty).get
  val prepareCall: TaskCall = wdl.workflow.findCallByName("do_prepare").get.asInstanceOf[TaskCall]
  val scatterCall: TaskCall = wdl.workflow.findCallByName("do_scatter").get.asInstanceOf[TaskCall]
  
  def makeKey(call: TaskCall, executionStatus: ExecutionStatus)(index: Int) = {
    BackendJobDescriptorKey(call, Option(index), 1) -> executionStatus
  }
  
  // Generates numbers from 1000 to 10000 with 1000 gap:
  // 1000, 2000, ..., 10000
  val sizes: Gen[Int] = Gen.range("size")(1000, 10000, 1000)
  
  // Generates executionStores using the given above sizes
  // Each execution store contains X simulated shards of "prepareCall" in status Done and X simulated shards of "scatterCall" in status NotStarted
  // This provides a good starting point to evaluate the speed of "runnableCalls", as it needs to iterate over all "NotStarted" keys, and for each one
  // look for their upstreams keys in status "Done"
  val executionStores: Gen[ExecutionStore] = for {
    size <- sizes
    doneMap = (0 until size map makeKey(prepareCall, ExecutionStatus.Done)).toMap
    notStartedMap = (0 until size map makeKey(scatterCall, ExecutionStatus.NotStarted)).toMap
    finalMap: Map[JobKey, ExecutionStatus] = doneMap ++ notStartedMap
  } yield new ExecutionStore(finalMap, true)
  
  performance of "ExecutionStore" in {
    measure method "runnableCalls" in {
      using(executionStores) in { es =>
        es.runnableScopes
      }
    }
  }
}
