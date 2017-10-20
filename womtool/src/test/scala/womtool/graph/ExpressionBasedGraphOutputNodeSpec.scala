package womtool.graph

import cats.data.Validated.{Invalid, Valid}
import wom.types.WomIntegerType
import wom.expression.PlaceholderWomExpression
import wom.graph._

class ExpressionBasedGraphOutputNodeSpec extends WomDotGraphTest {

  val expressionOutputGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType)
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType)

    // Declare a port output from i:
    val jOutput = PortBasedGraphOutputNode(WomIdentifier("j_out"), WomIntegerType, jInputNode.singleOutputPort)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression output using both i and j:
    val xOutputValidation = ExpressionBasedGraphOutputNode.fromInputMapping(
      WomIdentifier("x_out"), ijExpression, ijExpression.fixedWomType, Map(
      "i" -> iInputNode.singleOutputPort,
      "j" -> jInputNode.singleOutputPort))

    import common.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = xOutputValidation flatMap { xOutput => Graph.validateAndConstruct(Set(iInputNode, jInputNode, jOutput, xOutput)) }

    graph match {
      case Valid(g) => g
      case Invalid(errors) => throw new Exception(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }
  }

  val expressionOutputDot =
    """digraph "ExpressionBasedGraphOutputNodes"
      |{
      |compound=true;
      |"PORT0" -> "PORT1"
      |  "PORT2" -> "PORT3"
      |  "PORT0" -> "PORT5"
      |
      |subgraph cluster_0 {
      |  style=filled;
      |  fillcolor=lightskyblue1;
      |  "NODE7" [shape=plaintext label="Int i"]
      |  "PORT2" [shape=hexagon label="i"];
      |}
      |
      |
      |subgraph cluster_1 {
      |  style=filled;
      |  fillcolor=lightskyblue1;
      |  "NODE9" [shape=plaintext label="Int j"]
      |  "PORT0" [shape=hexagon label="j"];
      |}
      |
      |
      |subgraph cluster_2 {
      |  style=filled;
      |  fillcolor=palegreen;
      |  "NODE11" [shape=plaintext label="Int j_out"]
      |  "PORT1" [shape=oval label="j_out"];
      |}
      |
      |
      |subgraph cluster_3 {
      |  style=filled;
      |  fillcolor=palegreen;
      |  "NODE14" [shape=plaintext label="Int x_out"]
      |  "PORT3" [shape=oval label="i"];
      |  "PORT5" [shape=oval label="j"];
      |}
      |
      |}
      |""".stripMargin

  // TODO WOM uncomment
  // override val cases = List(WomDotGraphTestCase("ExpressionBasedGraphOutputNodes", expressionOutputGraph, expressionOutputDot))
  override val cases = List.empty

  tests()
}
