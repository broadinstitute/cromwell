package wdl4s.wdl.wom

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.wom.WdlWomExpressionsAsInputsSpec.Wdl
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.graph.TaskCallNode

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

    val workflowGraph = namespace.womExecutable.flatMap(_.graph) match {
      case Valid(g) => g
      case Invalid(errors) => fail(s"Unable to build wom graph from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    val callNodes = workflowGraph.nodes.toList collect { case callNode: TaskCallNode => callNode.name -> callNode } toMap

    val c = callNodes("c")
    c.expressionBasedInputs should have size 1
    val inputExpression = c.expressionBasedInputs.values.head

    List("a", "b") foreach { x =>
      val connectedInputPort = inputExpression.inputMapping(s"$x.int_out")
      val upstream = connectedInputPort.upstream
      upstream.name shouldBe "int_out"
      val upstreamTaskCallNode = upstream.graphNode.asInstanceOf[TaskCallNode]
      upstreamTaskCallNode.name shouldBe x
      // Instance equality test
      (upstreamTaskCallNode eq callNodes(x)) shouldBe true
    }
  }
}
