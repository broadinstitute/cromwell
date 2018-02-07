package womtool.graph

import wdl.WdlNamespaceWithWorkflow
import wdl.transforms.draft2.wdlom2wom._
import wom.transforms.WomExecutableMaker.ops._

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
    namespace.toWomExecutable() match {
      case Right(executable) => executable.graph
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
      |    "NODE2" [shape=plaintext label="call tsk"]
      |    "PORT0" [shape=hexagon label="Array[Array[File]] out"];
      |    "PORT3" [shape=hexagon label="Array[File] reads_1"];
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
