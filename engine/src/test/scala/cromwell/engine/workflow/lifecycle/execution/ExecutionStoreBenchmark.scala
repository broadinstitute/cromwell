package cromwell.engine.workflow.lifecycle.execution

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{apply => _, _}
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.keys.ScatterCollectorKey
import cromwell.engine.workflow.lifecycle.execution.stores.ActiveExecutionStore
import cromwell.util.SampleWdl
import org.scalameter.api._
import org.scalameter.picklers.Implicits._
import spray.json.DefaultJsonProtocol
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wom.graph.{CommandCallNode, ScatterNode}
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wom.expression.NoIoFunctionSet
import wom.transforms.WomExecutableMaker.ops._

/**
  * Benchmarks the performance of the execution store using ScalaMeter (http://scalameter.github.io/)
  * This is not run automatically by "sbt test". To run this test specifically, either use intellij integration, or run
  * sbt "engine/benchmark:testOnly cromwell.engine.workflow.lifecycle.execution.ExecutionStoreBenchmark"
  * sbt benchmark:test will run all ScalaMeter tests
  */
object ExecutionStoreBenchmark extends Bench[Double] with DefaultJsonProtocol {

  import spray.json._
  /* Benchmark configuration */
  lazy val measurer = new Measurer.Default
  lazy val executor = LocalExecutor(new Executor.Warmer.Default, Aggregator.average, measurer)
  lazy val reporter = new LoggingReporter[Double]
  lazy val persistor = Persistor.None
  
  val inputJson = Option(SampleWdl.PrepareScatterGatherWdl().rawInputs.toJson.compactPrint)
  val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.PrepareScatterGatherWdl().workflowSource(), Seq.empty).get
  val graph = namespace.toWomExecutable(inputJson, NoIoFunctionSet, strictValidation = true).getOrElse(throw new Exception("Failed to build womExecutable")).graph
  val prepareCall: CommandCallNode = graph.calls.find(_.localName == "do_prepare").get.asInstanceOf[CommandCallNode]
  val scatterCall: CommandCallNode = graph.allNodes.find(_.localName == "do_scatter").get.asInstanceOf[CommandCallNode]
  val scatter: ScatterNode = graph.scatters.head
  
  private def makeKey(call: CommandCallNode, executionStatus: ExecutionStatus)(index: Int) = {
    BackendJobDescriptorKey(call, Option(index), 1) -> executionStatus
  }

  // Generates executionStores using the given above sizes
  // Each execution store contains X simulated shards of "prepareCall" in status Done and X simulated shards of "scatterCall" in status NotStarted
  // This provides a good starting point to evaluate the speed of "runnableCalls", as it needs to iterate over all "NotStarted" keys, and for each one
  // look for their upstreams keys in status "Done"
  private def stores(sizes: Gen[Int]): Gen[ActiveExecutionStore] = {
    for {
      size <- sizes
      doneMap = (0 until size map makeKey(prepareCall, ExecutionStatus.Done)).toMap
      collectorKeys = scatter.outputMapping.map(om => ScatterCollectorKey(om, size, ScatterNode.DefaultScatterCollectionFunction) -> ExecutionStatus.NotStarted).toMap
      notStartedMap = (0 until size map makeKey(scatterCall, ExecutionStatus.NotStarted)).toMap ++ collectorKeys
      finalMap: Map[JobKey, ExecutionStatus] = doneMap ++ notStartedMap
    } yield ActiveExecutionStore(finalMap, needsUpdate = true)
  }
  
  performance of "ExecutionStore" in {
    // Measures how fast the execution store can find runnable calls with lots of "Done" calls and "NotStarted" calls.
    // Other "shapes" would be valuable to get a better sense of how this method behaves in various situations (with Collector Keys etc...)  
    measure method "update" in {
      val sizes: Gen[Int] = Gen.range("size")(from = 1000, upto = 10000, hop = 1000)
      using(stores(sizes)) in { es =>
        es.update
      }
    }

    measure method "backendJobDescriptorKeyForNode" in {
      val sizes: Gen[Int] = Gen.range("size")(from = 100000, upto = 1000000, hop = 100000)
      using(stores(sizes)) in { es =>
        es.backendJobDescriptorKeyForNode(scatterCall)
      }
    }
  }
}
