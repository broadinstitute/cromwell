package wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.expression._
import wom.types.{WdlIntegerType, WdlStringType}

class GraphOutputNodeSpec extends FlatSpec with Matchers {

  behavior of "ExpressionBasedGraphOutputNode"

  it should "construct an ExpressionBasedGraphOutputNode node if inputs are available" in {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WdlIntegerType)
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WdlIntegerType)

    // Declare a port output from i:
    val jOutput = PortBasedGraphOutputNode(WomIdentifier("j_out"), WdlIntegerType, jInputNode.singleOutputPort)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    // Declare the expression output using both i and j:
    val xOutputValidation = ExpressionBasedGraphOutputNode.linkWithInputs(WomIdentifier("x_out"), WdlStringType, ijExpression, Map(
      "i" -> iInputNode.singleOutputPort,
      "j" -> jInputNode.singleOutputPort))

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = xOutputValidation flatMap { xOutput => Graph.validateAndConstruct(Set(iInputNode, jInputNode, jOutput, xOutput)) }

    graph match {
      case Valid(g) => validate(g)
      case Invalid(errors) => fail(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }

    def validate(graph: Graph) = {
      val expressionOutput = graph.nodes.find {
        case g: GraphOutputNode => g.localName == "x_out"
        case _ => false
      }

      expressionOutput.isDefined should be(true)
      expressionOutput.get.upstream should be(Set(iInputNode, jInputNode))

      val portOutput = graph.nodes.find {
        case g: GraphOutputNode => g.localName == "j_out"
        case _ => false
      }

      portOutput.isDefined should be(true)
      portOutput.get.upstream should be(Set(jInputNode))
    }
  }

}
