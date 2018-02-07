package wdl.transforms

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow, WdlWorkflow}
import wom.graph.CommandCallNode
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExpressionNode
import WdlWomExpressionsAsInputsSpec._
import wom.transforms.WomWorkflowDefinitionMaker
import wom.transforms.WomWorkflowDefinitionMaker.ops._

import scala.language.postfixOps

object WdlWomExpressionsAsInputsSpec {
  val Wdl =
    // Using calls as inputs since declaration inputs are currently not supported.
    """
      |workflow foo {
      |    call x as a
      |    call x as b
      |
      |    call x as c { input: int_in = a.int_out + b.int_out }
      |}
      |
      |task x {
      |    Int int_in
      |    command {}
      |    output {
      |        Int int_out = int_in
      |    }
      |}
    """.stripMargin
}


class WdlWomExpressionsAsInputsSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {
  behavior of "WdlWomExpressionsAsInputs"

  it should "wire up input expressions for a WDL workflow" in {

    val namespace = WdlNamespace.loadUsingSource(Wdl, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    val workflowGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph) match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom graph from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    val callNodes = workflowGraph.nodes.toList collect { case callNode: CommandCallNode => callNode.localName -> callNode } toMap

    val a = callNodes("a")
    val b = callNodes("b")
    val c = callNodes("c")
    c.inputPorts should have size 1
    val cInputExpressionNode = c.inputPorts.map(_.upstream).head.graphNode.asInstanceOf[ExpressionNode]
    cInputExpressionNode.inputPorts.map(_.upstream) shouldBe a.outputPorts ++ b.outputPorts

    val inputExpression = c.inputDefinitionMappings
      .head._2.select[OutputPort].get
      .graphNode.asInstanceOf[ExpressionNode]  

    List("a", "b") foreach { x =>
      val connectedInputPort = inputExpression.inputMapping(s"$x.int_out")
      val upstream = connectedInputPort.upstream
      upstream.name shouldBe "int_out"
      val upstreamTaskCallNode = upstream.graphNode.asInstanceOf[CommandCallNode]
      upstreamTaskCallNode.localName shouldBe x
      // Instance equality test
      (upstreamTaskCallNode eq callNodes(x)) shouldBe true
    }
  }
}
