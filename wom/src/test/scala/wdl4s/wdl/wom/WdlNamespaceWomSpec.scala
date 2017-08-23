package wdl4s.wdl.wom

import lenthall.collections.EnhancedCollections._
import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.graph.{CallNode, GraphInputNode, PortBasedGraphOutputNode, TaskCallNode}

class WdlNamespaceWomSpec extends FlatSpec with Matchers {
  
  "A WdlNamespace for 3step" should "provide conversion to WOM" in {
    val threeStep =
      """
        |task ps {
        |  command {
        |    ps
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |}
        |
        |task cgrep {
        |  String pattern
        |  File in_file
        |
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |task wc {
        |  File in_file
        |  command {
        |    cat ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |workflow three_step {
        |  call ps
        |  call cgrep {
        |    input: in_file = ps.procs
        |  }
        |  call wc {
        |    input: in_file = ps.procs
        |  }
        |}
      """.stripMargin

    val namespace = WdlNamespace.loadUsingSource(threeStep, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val wom3Step = namespace.womExecutable.flatMap(_.graph)
    
    val workflowGraph = wom3Step match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom version of 3step from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }
    
    workflowGraph.nodes collect { case gin: GraphInputNode => gin.name } should be(Set("cgrep.pattern"))
    workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon.name } should be(Set("wc.count", "cgrep.count", "ps.procs"))
    
    val ps: TaskCallNode = workflowGraph.nodes.collectFirst({ case ps: TaskCallNode if ps.name == "ps" => ps }).get
    val cgrep: TaskCallNode = workflowGraph.nodes.collectFirst({ case cgrep: TaskCallNode if cgrep.name == "cgrep" => cgrep }).get
    val cgrepPatternInput = workflowGraph.nodes.collectFirst({ case cgrepInput: GraphInputNode if cgrepInput.name == "cgrep.pattern" => cgrepInput }).get
    val wc: TaskCallNode = workflowGraph.nodes.collectFirst({ case wc: TaskCallNode if wc.name == "wc" => wc }).get

    workflowGraph.nodes.filterByType[CallNode] should be(Set(ps, cgrep, wc))
    ps.inputPorts.map(_.name) should be(Set.empty)
    cgrep.inputPorts.map(_.name) should be(Set("pattern", "ps.procs"))
    wc.inputPorts.map(_.name) should be(Set("ps.procs"))

    ps.upstream shouldBe empty
    cgrep.upstream shouldBe Set(ps, cgrepPatternInput)
    wc.upstream shouldBe Set(ps)
  }

}
