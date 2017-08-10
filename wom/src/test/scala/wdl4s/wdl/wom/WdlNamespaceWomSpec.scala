package wdl4s.wdl.wom

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.graph.{CallNode, GraphInputNode, PortBasedGraphOutputNode}

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
    val wom3Step = namespace.womExecutable
    
    val workflowGraph = wom3Step.graph match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom version of 3step from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }
    
    workflowGraph.nodes collect { case gin: GraphInputNode => gin.name } should be(Set("cgrep.pattern"))
    workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon.name } should be(Set("wc.count", "cgrep.count", "ps.procs"))
    workflowGraph.nodes collect { case cn: CallNode => cn.name } should be(Set("wc", "cgrep", "ps"))
    
    val ps = workflowGraph.nodes.collectFirst({ case ps: CallNode if ps.name == "ps" => ps }).get
    val cgrep = workflowGraph.nodes.collectFirst({ case cgrep: CallNode if cgrep.name == "cgrep" => cgrep }).get
    val cgrepPatternInput = workflowGraph.nodes.collectFirst({ case cgrepInput: GraphInputNode if cgrepInput.name == "cgrep.pattern" => cgrepInput }).get
    val wc = workflowGraph.nodes.collectFirst({ case wc: CallNode if wc.name == "wc" => wc }).get
    
    ps.upstream shouldBe empty
    cgrep.upstream shouldBe Set(ps, cgrepPatternInput)
    wc.upstream shouldBe Set(ps)
  }

}
