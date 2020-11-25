package wom.graph

import cats.data.Validated.{Invalid, Valid}
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.expression._
import wom.types.{WomIntegerType, WomStringType}


class GraphOutputNodeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "ExpressionBasedGraphOutputNode"

  it should "construct an ExpressionBasedGraphOutputNode node if inputs are available" in {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType, "i")
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType, "j")

    // Declare a port output from i:
    val jOutput = PortBasedGraphOutputNode(WomIdentifier("j_out"), WomIntegerType, jInputNode.singleOutputPort)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression output using both i and j:
    val xOutputValidation = ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier("x_out"), ijExpression, WomStringType, Map(
      "i" -> iInputNode.singleOutputPort,
      "j" -> jInputNode.singleOutputPort))

    import common.validation.ErrorOr.ShortCircuitingFlatMap
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
