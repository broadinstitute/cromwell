package wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.expression._
import wom.types.WdlIntegerType

class ExpressionNodeSpec extends FlatSpec with Matchers {

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
  it should "construct an ExpressionNode if inputs are available" in {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WdlIntegerType)
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WdlIntegerType)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    // Declare the expression node using both i and j:
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      xDeclarationNode <- ExpressionNode.linkWithInputs(WomIdentifier("x"), ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))
      xOutputNode = PortBasedGraphOutputNode(WomIdentifier("x_out"), WdlIntegerType, xDeclarationNode.singleExpressionOutputPort)
      g <- Graph.validateAndConstruct(Set(iInputNode, jInputNode, xDeclarationNode, xOutputNode))
    } yield g

    graph match {
      case Valid(g) => validate(g)
      case Invalid(errors) => fail(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }

    def validate(graph: Graph) = {
      val x_outGraphOutputNode = graph.nodes.find {
        case g: GraphOutputNode => g.localName == "x_out"
        case _ => false
      }.getOrElse(fail("No 'x_out' GraphOutputNode in the graph"))

      x_outGraphOutputNode.upstream.size should be(1)
      val xExpressionNode = x_outGraphOutputNode.upstream.head.asInstanceOf[ExpressionNode]
      xExpressionNode.localName should be("x")
      xExpressionNode.womType should be(WdlIntegerType)

      xExpressionNode.upstream should be(Set(iInputNode, jInputNode))
    }
  }

}
