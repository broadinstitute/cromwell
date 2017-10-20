package wdltool.graph

import cats.data.Validated.{Invalid, Valid}
import wdl4s.wdl.types.WdlIntegerType
import wdl4s.wom.expression.PlaceholderWomExpression
import wdl4s.wom.graph._

class ExpressionNodeSpec extends WomDotGraphTest {

  val expressionNodeGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode("i", WdlIntegerType)
    val jInputNode = RequiredGraphInputNode("j", WdlIntegerType)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WdlIntegerType)

    // Declare the expression node using both i and j:
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      xDeclarationNode <- ExpressionNode.linkWithInputs("x", ijExpression, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))
      xOutputNode = PortBasedGraphOutputNode("x_out", WdlIntegerType, xDeclarationNode.singleExpressionOutputPort)
      g <- Graph.validateAndConstruct(Set(iInputNode, jInputNode, xDeclarationNode, xOutputNode))
    } yield g

    graph match {
      case Valid(g) => g
      case Invalid(errors) => throw new Exception(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }
  }

  val expressionNodeDot =
    """digraph "ExpressionNodes"
      |{
      |compound=true;
      |"PORT0" -> "PORT1"
      |"PORT2" -> "PORT3"
      |"PORT4" -> "PORT5"
      |subgraph cluster_0 {
      |style=filled;
      |fillcolor=lightskyblue1;
      |"NODE6" [shape=plaintext label="Int i"]
      |"PORT0" [shape=hexagon label="i"];
      |}
      |subgraph cluster_1 {
      |style=filled;
      |fillcolor=lightskyblue1;
      |"NODE7" [shape=plaintext label="Int j"]
      |"PORT2" [shape=hexagon label="j"];
      |}
      |subgraph cluster_2 {
      |style=filled;
      |fillcolor=white;
      |"NODE8" [shape=plaintext label="Int x = f(i, j)"]
      |"PORT4" [shape=hexagon label="x"];
      |"PORT1" [shape=oval label="i"];
      |"PORT3" [shape=oval label="j"];
      |}
      |subgraph cluster_3 {
      |style=filled;
      |fillcolor=palegreen;
      |"NODE9" [shape=plaintext label="Int x_out"]
      |"PORT5" [shape=oval label="x_out"];
      |}
      |}
      |""".stripMargin

  override val cases = List(WomDotGraphTestCase("ExpressionNodes", expressionNodeGraph, expressionNodeDot))

  tests()
}
