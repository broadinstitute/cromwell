package womtool.graph

import cats.data.Validated.{Invalid, Valid}
import wom.expression.PlaceholderWomExpression
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomIntegerType

class ExpressionNodeSpec extends WomDotGraphTest {

  val expressionNodeGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType)
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression node using both i and j:
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      xDeclarationNode <- ExposedExpressionNode.fromInputMapping(
        WomIdentifier("x"), ijExpression, WomIntegerType, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))
      xOutputNode = PortBasedGraphOutputNode(WomIdentifier("x_out"), WomIntegerType, xDeclarationNode.singleExpressionOutputPort)
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

  // TODO WOM uncomment
  // override val cases = List(WomDotGraphTestCase("ExpressionNodes", expressionNodeGraph, expressionNodeDot))
  override val cases = List.empty

  tests()
}
