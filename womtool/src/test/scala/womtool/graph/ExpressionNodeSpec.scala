package womtool.graph

import cats.data.Validated.{Invalid, Valid}
import wom.expression.PlaceholderWomExpression
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomIntegerType

class ExpressionNodeSpec extends WomDotGraphTest {

  val expressionNodeGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType, "i")
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType, "j")

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression node using both i and j:
    import common.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = for {
      xDeclarationNode <- ExposedExpressionNode.fromInputMapping(
        WomIdentifier("x"), ijExpression, WomIntegerType, Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort))
      xOutputNode = PortBasedGraphOutputNode(WomIdentifier("x_out"), WomIntegerType, xDeclarationNode.singleOutputPort)
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
      |  compound=true;
      |  "PORT0" -> "PORT1"
      |  "PORT2" -> "PORT3"
      |  "PORT4" -> "PORT5"
      |  subgraph cluster_0 {
      |    style="filled,solid";
      |    fillcolor=lightskyblue1;
      |    "NODE6" [shape=plaintext label="Int i"]
      |    "PORT0" [shape=hexagon label="Int i"];
      |  }
      |  subgraph cluster_1 {
      |    style="filled,solid";
      |    fillcolor=lightskyblue1;
      |    "NODE7" [shape=plaintext label="Int j"]
      |    "PORT2" [shape=hexagon label="Int j"];
      |  }
      |  subgraph cluster_2 {
      |    style="filled,solid";
      |    fillcolor=white;
      |    "NODE8" [shape=plaintext label="Int x = f(i, j)"]
      |    "PORT4" [shape=hexagon label="Int x"];
      |    "PORT1" [shape=oval label="Int i"];
      |    "PORT3" [shape=oval label="Int j"];
      |  }
      |  subgraph cluster_3 {
      |    style="filled,solid";
      |    fillcolor=yellowgreen;
      |    "NODE9" [shape=plaintext label="Int x_out"]
      |    "PORT5" [shape=oval label="Int x_out"];
      |  }
      |}
      |""".stripMargin

  override val cases = List(WomDotGraphTestCase("ExpressionNodes", expressionNodeGraph, expressionNodeDot))

  tests()
}
