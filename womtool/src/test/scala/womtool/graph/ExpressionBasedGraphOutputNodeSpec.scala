package womtool.graph

import cats.data.Validated.{Invalid, Valid}
import wom.types.WomIntegerType
import wom.expression.PlaceholderWomExpression
import wom.graph._

class ExpressionBasedGraphOutputNodeSpec extends WomDotGraphTest {

  val expressionOutputGraph = {
    // Two inputs:
    val iInputNode = RequiredGraphInputNode(WomIdentifier("i"), WomIntegerType, "i")
    val jInputNode = RequiredGraphInputNode(WomIdentifier("j"), WomIntegerType, "j")

    // Declare a port output from i:
    val jOutput = PortBasedGraphOutputNode(WomIdentifier("j_out"), WomIntegerType, jInputNode.singleOutputPort)

    // Declare an expression that needs both an "i" and a "j":
    val ijExpression = PlaceholderWomExpression(Set("i", "j"), WomIntegerType)

    // Declare the expression output using both i and j:
    val xOutputValidation = ExpressionBasedGraphOutputNode.fromInputMapping(
      WomIdentifier("x_out"),
      ijExpression,
      ijExpression.fixedWomType,
      Map("i" -> iInputNode.singleOutputPort, "j" -> jInputNode.singleOutputPort)
    )

    import common.validation.ErrorOr.ShortCircuitingFlatMap
    val graph = xOutputValidation flatMap { xOutput =>
      Graph.validateAndConstruct(Set(iInputNode, jInputNode, jOutput, xOutput))
    }

    graph match {
      case Valid(g) => g
      case Invalid(errors) => throw new Exception(s"Unable to build WOM graph: ${errors.toList.mkString(" ")}")
    }
  }

  val expressionOutputDot =
    """digraph "ExpressionBasedGraphOutputNodes"
      |{
      |  compound=true;
      |  "PORT0" -> "PORT1"
      |  "PORT2" -> "PORT3"
      |  "PORT0" -> "PORT4"
      |  subgraph cluster_0 {
      |    style="filled,solid";
      |    fillcolor=lightskyblue1;
      |    "NODE5" [shape=plaintext label="Int i"]
      |    "PORT2" [shape=hexagon label="Int i"];
      |  }
      |  subgraph cluster_1 {
      |    style="filled,solid";
      |    fillcolor=lightskyblue1;
      |    "NODE6" [shape=plaintext label="Int j"]
      |    "PORT0" [shape=hexagon label="Int j"];
      |  }
      |  subgraph cluster_2 {
      |    style="filled,solid";
      |    fillcolor=yellowgreen;
      |    "NODE7" [shape=plaintext label="Int j_out"]
      |    "PORT1" [shape=oval label="Int j_out"];
      |  }
      |  subgraph cluster_3 {
      |    style="filled,solid";
      |    fillcolor=palegreen;
      |    "NODE8" [shape=plaintext label="Int x_out"]
      |    "PORT9" [shape=hexagon label="Int x_out"];
      |    "PORT3" [shape=oval label="Int i"];
      |    "PORT4" [shape=oval label="Int j"];
      |  }
      |}
      |""".stripMargin

  override val cases = List(
    WomDotGraphTestCase("ExpressionBasedGraphOutputNodes", expressionOutputGraph, expressionOutputDot)
  )
  tests()
}
