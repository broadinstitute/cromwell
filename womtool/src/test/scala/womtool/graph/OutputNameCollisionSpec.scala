package womtool.graph

import common.collections.EnhancedCollections._
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers._
import wom.callable.WorkflowDefinition
import wom.transforms.WomBundleMaker.ops._

class OutputNameCollisionSpec extends WomDotGraphTest {

  val wdl =
    """
      |workflow wf {
      |  # This call has an output with a reference to 'out'. That **must not** be resolved to the wf.out!
      |  call tsk
      |  output {
      |    Array[File] out = tsk.out[0]
      |  }
      |}
      |
      |task tsk {
      |  command {
      |    # ...
      |  }
      |  output {
      |    Array[Array[File]] out = [["out.txt"]]
      |    Array[File] reads_1 = out[0]
      |  }
      |}
    """.stripMargin

  val outputCollisionWdlGraph = {
    val namespace = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

    namespace.toWomBundle match {
      case Right(bundle) => (bundle.allCallables.filterByType[WorkflowDefinition]: Set[WorkflowDefinition]).head.graph
      case Left(errors) => throw new Exception(errors.toList.mkString(", "))
    }
  }

  val outputCollisionWdlDot =
    """digraph "non_colliding_output_names"
      |{
      |  compound=true;
      |  "PORT0" -> "PORT1"
      |
      |  subgraph cluster_0 {
      |    style="filled,solid";
      |    fillcolor=white;
      |    "NODE2" [shape=plaintext label="call wf.tsk (tsk)"]
      |    "PORT0" [shape=hexagon label="Array[Array[File]] tsk.out"];
      |    "PORT3" [shape=hexagon label="Array[File] tsk.reads_1"];
      |  }
      |  subgraph cluster_1 {
      |    style="filled,solid";
      |    fillcolor=palegreen;
      |    "NODE4" [shape=plaintext label="Array[File] out"]
      |    "PORT5" [shape=hexagon label="Array[File] out"];
      |    "PORT1" [shape=oval label="Array[Array[File]] tsk.out"];
      |  }
      |}
      |""".stripMargin

  override val cases = List(WomDotGraphTestCase("non_colliding_output_names", outputCollisionWdlGraph, outputCollisionWdlDot))

  tests()
}
