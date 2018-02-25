package wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.expression._
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode, PlainAnonymousExpressionNode}
import wom.types.WomIntegerType

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
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType, "i")
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType, "j")

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression node using both i and j:
    import common.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      xDeclarationNode <- AnonymousExpressionNode.fromInputMapping(WomIdentifier("x"), ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort), PlainAnonymousExpressionNode.apply)
      xOutputNode = PortBasedGraphOutputNode(WomIdentifier("x_out"), WomIntegerType, xDeclarationNode.singleOutputPort)
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
      xExpressionNode.womType should be(WomIntegerType)

      xExpressionNode.upstream should be(Set(iInputNode, jInputNode))
    }
  }

}
