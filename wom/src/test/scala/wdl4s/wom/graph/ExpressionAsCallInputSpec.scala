package wdl4s.wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.types.WdlIntegerType
import wdl4s.wom.callable.TaskDefinitionSpec
import wdl4s.wom.expression._
import wdl4s.wom.graph.CallNode.CallWithInputs

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
    val iInputNode = RequiredGraphInputNode("i", WdlIntegerType)
    val jInputNode = RequiredGraphInputNode("j", WdlIntegerType)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    def validateCallResult(callWithInputs: CallWithInputs) = {
      callWithInputs.inputs should be(Set.empty)
      callWithInputs.call.upstream should be(Set(iInputNode, jInputNode))
    }

    // Use that as an input to a one-input task:
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      callNode <- CallNode.callWithInputs("foo", TaskDefinitionSpec.oneInputTask, Map.empty, Set(GraphNodeInputExpression("bar", ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))))
      _ = validateCallResult(callNode)
      g <- Graph.validateAndConstruct(Set(iInputNode, jInputNode, callNode.call))
    } yield g

    graph match {
      case Valid(_) => // Great!
      case Invalid(errors) => fail(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }
  }

}
