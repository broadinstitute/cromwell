package wom

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.WdlWomExpressionsAsInputsSpec.Wdl
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{ExpressionNode, TaskCallNode}

import scala.language.postfixOps


object WdlWomExpressionsAsInputsSpec {
  val Wdl =
    // Using calls as inputs since declaration inputs are currently not supported.
    """
      |workflow foo {
      |    call a
      |    call b
      |
      |    call c { input: int_in = a.int_out + b.int_out }
      |}
      |
      |task a {
      |    Int int_in
      |    command {}
      |    output {
      |        Int int_out = int_in
      |    }
      |}
      |
      |task b {
      |    Int int_in
      |    command {}
      |    output {
      |        Int int_out = int_in
      |    }
      |}
      |
      |task c {
      |    Int int_in
      |    command {}
      |    output {
      |        Int int_out = int_in
      |    }
      |}
    """.stripMargin
}


class WdlWomExpressionsAsInputsSpec extends FlatSpec with Matchers {
  behavior of "WdlWomExpressionsAsInputs"

  it should "wire up input expressions for a WDL workflow" in {

    val namespace = WdlNamespace.loadUsingSource(Wdl, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val workflowGraph = namespace.workflow.womDefinition.flatMap(_.graph) match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom graph from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    val callNodes = workflowGraph.nodes.toList collect { case callNode: TaskCallNode => callNode.localName -> callNode } toMap

    val a = callNodes("a")
    val b = callNodes("b")
    val c = callNodes("c")
    c.inputPorts should have size 1
    val cInputExpressionNode = c.inputPorts.map(_.upstream).head.graphNode.asInstanceOf[ExpressionNode]
    cInputExpressionNode.inputPorts.map(_.upstream) shouldBe a.outputPorts ++ b.outputPorts

    val inputExpression = c.inputDefinitionMappings
      .head._2.select[OutputPort].get
      .graphNode.asInstanceOf[ExpressionNode].instantiatedExpression  

    List("a", "b") foreach { x =>
      val connectedInputPort = inputExpression.inputMapping(s"$x.int_out")
      val upstream = connectedInputPort.upstream
      upstream.name shouldBe "int_out"
      val upstreamTaskCallNode = upstream.graphNode.asInstanceOf[TaskCallNode]
      upstreamTaskCallNode.localName shouldBe x
      // Instance equality test
      (upstreamTaskCallNode eq callNodes(x)) shouldBe true
    }
  }
}
