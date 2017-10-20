package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{apply => _, _}
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.keys.ScatterCollectorKey
import cromwell.util.SampleWdl
import org.scalameter.api._
import org.scalameter.picklers.Implicits._
import spray.json.DefaultJsonProtocol
import wdl.WdlNamespaceWithWorkflow
import wom.graph.{ScatterNode, TaskCallNode}

/**
  * Benchmarks the performance of the execution store using ScalaMeter (http://scalameter.github.io/)
  * This is not run automatically by "sbt test". To run this test specifically, either use intellij integration, or run
  * sbt "project engine"  "benchmark:test-only cromwell.engine.workflow.lifecycle.execution.ExecutionStoreBenchmark"
  * sbt benchmark:test will run all ScalaMeter tests
  */
object ExecutionStoreBenchmark extends Bench[Double] with DefaultJsonProtocol {

  import spray.json._
  /* Benchmark configuration */
  lazy val measurer = new Measurer.Default
  lazy val executor = SeparateJvmsExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  lazy val reporter = new LoggingReporter[Double]
  lazy val persistor = Persistor.None
  
  val inputJson = Option(SampleWdl.PrepareScatterGatherWdl().rawInputs.toJson.compactPrint)
  val wdl = WdlNamespaceWithWorkflow.load(SampleWdl.PrepareScatterGatherWdl().workflowSource(), Seq.empty).get
  val graph = wdl.womExecutable(inputJson).getOrElse(throw new Exception("Failed to build womExecutable"))
    .graph.getOrElse(throw new Exception("Failed to build wom graph"))
  val prepareCall: TaskCallNode = graph.calls.find(_.localName == "do_prepare").get.asInstanceOf[TaskCallNode]
  val scatterCall: TaskCallNode = graph.allNodes.find(_.localName == "do_scatter").get.asInstanceOf[TaskCallNode]
  val scatter: ScatterNode = graph.scatters.head
  
  private def makeKey(call: TaskCallNode, executionStatus: ExecutionStatus)(index: Int) = {
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
    collectorKey = Map(ScatterCollectorKey(scatterCall, scatter.outputMapping, scatter, size) -> ExecutionStatus.NotStarted)
    notStartedMap = (0 until size map makeKey(scatterCall, ExecutionStatus.NotStarted)).toMap ++ collectorKey
    finalMap: Map[JobKey, ExecutionStatus] = doneMap ++ notStartedMap
  } yield new ExecutionStore(finalMap, true)
  
  performance of "ExecutionStore" in {
    // Measures how fast the execution store can find runnable calls with lots of "Done" calls and "NotStarted" calls.
    // Other "shapes" would be valuable to get a better sense of how this method behaves in various situations (with Collector Keys etc...)  
    measure method "runnableCalls" in {
      using(executionStores) in { es =>
        es.runnableScopes
      }
    }
  }
}
