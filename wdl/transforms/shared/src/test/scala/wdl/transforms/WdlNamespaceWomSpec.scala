package wdl.transforms

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow, WdlWomExpression, WdlWorkflow}
import wom.transforms.WomWorkflowDefinitionMaker
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExpressionNode
import common.collections.EnhancedCollections._
import wom.transforms.WomWorkflowDefinitionMaker.ops._

class WdlNamespaceWomSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {

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
    val wom3Step = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    val workflowGraph = wom3Step match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom version of 3step from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    val graphInputNodes = workflowGraph.nodes collect { case gin: ExternalGraphInputNode => gin }
    graphInputNodes should have size 1
    val patternInputNode = graphInputNodes.head
    patternInputNode.localName should be("cgrep.pattern")
    patternInputNode.fullyQualifiedName should be("three_step.cgrep.pattern")

    workflowGraph.nodes collect { case gon: ExpressionBasedGraphOutputNode => gon.localName } should be(Set("wc.count", "cgrep.count", "ps.procs"))

    val ps: CommandCallNode = workflowGraph.nodes.collectFirst({ case ps: CommandCallNode if ps.localName == "ps" => ps }).get
    val cgrep: CommandCallNode = workflowGraph.nodes.collectFirst({ case cgrep: CommandCallNode if cgrep.localName == "cgrep" => cgrep }).get
    val cgrepInFileExpression = {
      workflowGraph.nodes.collectFirst({ case cgrepInFile: ExpressionNode if cgrepInFile.localName == "cgrep.in_file" => cgrepInFile }).get
    }
    val wc: CommandCallNode = workflowGraph.nodes.collectFirst({ case wc: CommandCallNode if wc.localName == "wc" => wc }).get
    val wcInFileExpression = {
      workflowGraph.nodes.collectFirst({ case wcInFile: ExpressionNode if wcInFile.localName == "wc.in_file" => wcInFile }).get
    }

    workflowGraph.nodes.filterByType[CallNode] should be(Set(ps, cgrep, wc))
    ps.inputPorts.map(_.name) should be(Set.empty)
    cgrep.inputPorts.map(_.name) should be(Set("pattern", "in_file"))
    wc.inputPorts.map(_.name) should be(Set("in_file"))

    ps.upstream shouldBe empty
    cgrep.upstream shouldBe Set(cgrepInFileExpression, patternInputNode)
    wc.upstream shouldBe Set(wcInFileExpression)

    ps.inputDefinitionMappings shouldBe empty
    cgrep.inputDefinitionMappings should have size 2

    val cgrepFileInputDef = cgrep.callable.inputs.find(_.name == "in_file").get
    val cgrepInputs = cgrep.inputDefinitionMappings.toMap
    val inFileMapping = cgrepInputs(cgrepFileInputDef)
    inFileMapping.select[OutputPort].isDefined shouldBe true
    // This should be less ugly when we can access a string value from a womexpression
    inFileMapping.select[OutputPort].get
      .graphNode.asInstanceOf[ExpressionNode]
      .womExpression.asInstanceOf[WdlWomExpression]
      .wdlExpression.valueString shouldBe "ps.procs"

    val cgrepPatternInputDef = cgrep.callable.inputs.find(_.name == "pattern").get
    cgrepInputs(cgrepPatternInputDef).select[OutputPort].get eq patternInputNode.singleOutputPort shouldBe true
  }

}
