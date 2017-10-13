package wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.TaskDefinitionSpec
import wom.expression._
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold}
import wom.types.WdlIntegerType

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
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WdlIntegerType)
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WdlIntegerType)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    // Use that as an input to a one-input task:
    val expressionNode = ExpressionNode
      .linkWithInputs(WomIdentifier("bar"), ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))
      .getOrElse(fail("Failed to build expression node"))

    val callNodeBuilder = new CallNodeBuilder()

    val inputDefinition = TaskDefinitionSpec.oneInputTask.inputs.head
    
    val inputDefinitionFold = InputDefinitionFold(
      mappings = Map(inputDefinition -> expressionNode.inputDefinitionPointer),
      callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleExpressionOutputPort)),
      newExpressionNodes = Set(expressionNode)
    )
    
    val callNodeWithInputs = callNodeBuilder.build(
      WomIdentifier("foo"),
      TaskDefinitionSpec.oneInputTask,
      inputDefinitionFold
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
