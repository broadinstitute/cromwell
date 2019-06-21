package wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.CommandTaskDefinitionSpec
import wom.expression._
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold}
import wom.graph.expression.{AnonymousExpressionNode, TaskCallInputExpressionNode}
import wom.types.WomIntegerType

class ExpressionAsCallInputSpec extends FlatSpec with Matchers {

  behavior of "ExpressionBasedGraphOutputNode"

  /**
    * Roughly equivalent to the WDL (except that the expression could be anything):
    *
    * workflow foo {
    *   Int i
    *   Int j
    *   Int x = i + j
    *
    *   output {
    *     Int x_out = x
    *   }
    * }
    *
    */
  it should "create and wire in InstantiatedExpressions where appropriate" in {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType, "i")
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType, "j")

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Use that as an input to a one-input task:
    val expressionNode = AnonymousExpressionNode.fromInputMapping(
      WomIdentifier("bar"), ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort), TaskCallInputExpressionNode.apply)
      .getOrElse(fail("Failed to build expression node"))

    val callNodeBuilder = new CallNodeBuilder()

    val inputDefinition = CommandTaskDefinitionSpec.oneInputTask.inputs.head

    val inputDefinitionFold = InputDefinitionFold(
      mappings = List(inputDefinition -> expressionNode.inputDefinitionPointer),
      callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleOutputPort)),
      newExpressionNodes = Set(expressionNode)
    )

    val callNodeWithInputs = callNodeBuilder.build(
      WomIdentifier("foo"),
      CommandTaskDefinitionSpec.oneInputTask,
      inputDefinitionFold,
      Set.empty,
      None
    )

    def validateCallResult(callWithInputs: CallNodeAndNewNodes) = {
      callWithInputs.newInputs should be(Set.empty)
      callWithInputs.node.upstream should be(Set(expressionNode))
    }

    validateCallResult(callNodeWithInputs)

    val graph = Graph.validateAndConstruct(Set(iInputNode, jInputNode, callNodeWithInputs.node, expressionNode))

    graph match {
      case Valid(_) => // Great!
      case Invalid(errors) => fail(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }
  }

}
