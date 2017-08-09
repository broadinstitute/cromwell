package wdltool.graph

import cats.data.Validated.{Invalid, Valid}
import wdl4s.wdl.types.WdlIntegerType
import wdl4s.wom.expression.PlaceholderWomExpression
import wdl4s.wom.graph.{ExpressionBasedGraphOutputNode, Graph, PortBasedGraphOutputNode, RequiredGraphInputNode}

class ExpressionBasedGraphOutputNodeSpec extends WomDotGraphTest {

  val expressionOutputGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode("i", WdlIntegerType)
    val jInputNode = RequiredGraphInputNode("j", WdlIntegerType)

    // Declare a port output from i:
    val jOutput = PortBasedGraphOutputNode("j_out", WdlIntegerType, jInputNode.singleOutputPort)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    // Declare the expression output using both i and j:
    val xOutputValidation = ExpressionBasedGraphOutputNode.linkWithInputs("x_out", ijExpression, Map(
      "i" -> iInputNode.singleOutputPort,
      "j" -> jInputNode.singleOutputPort))

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
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
      |  "PORT2" [shape=hexagon label="i"];
      |  "NODE7" [shape=plaintext label="Int i"]
      |}
      |
      |
      |subgraph cluster_1 {
      |  style=filled;
      |  fillcolor=lightskyblue1;
      |  "PORT0" [shape=hexagon label="j"];
      |  "NODE9" [shape=plaintext label="Int j"]
      |}
      |
      |
      |subgraph cluster_2 {
      |  style=filled;
      |  fillcolor=palegreen;
      |  "PORT1" [shape=oval label="j_out"];
      |  "NODE11" [shape=plaintext label="Int j_out"]
      |}
      |
      |
      |subgraph cluster_3 {
      |  style=filled;
      |  fillcolor=palegreen;
      |  "PORT3" [shape=oval label="i"];
      |  "PORT5" [shape=oval label="j"];
      |  "NODE14" [shape=plaintext label="Int x_out"]
      |}
      |
      |}
      |""".stripMargin

  override val cases = List(WomDotGraphTestCase("ExpressionBasedGraphOutputNodes", expressionOutputGraph, expressionOutputDot))

  tests()
}
