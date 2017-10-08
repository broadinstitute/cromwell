package wdl

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import lenthall.validation.ErrorOr._
import wdl4s.parser.WdlParser.Ast
import wdl.types.WdlBooleanType
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph._
/**
  * Represents an If block in WDL
  *
  * @param index Index of the if block. The index is computed during tree generation to reflect WDL scope structure.
  * @param condition WDL Expression representing the condition in which to execute this If-block
  */
case class If(index: Int, condition: WdlExpression, ast: Ast) extends WdlGraphNodeWithUpstreamReferences with WorkflowScoped {
  val unqualifiedName = s"${If.FQNIdentifier}_$index"
  override def appearsInFqn = false

  final val upstreamReferences = condition.variableReferences

  override def toString: String = s"[If fqn=$fullyQualifiedName, condition=${condition.toWdlString}]"
}

object If {
  val FQNIdentifier = "$if"

  /**
    * @param index Index of the if block. The index is computed during tree generation to reflect WDL scope structure.
    */
  def apply(ast: Ast, index: Int): If = {
    new If(index, WdlExpression(ast.getAttribute("expression")), ast)
  }

  def womConditionalNode(ifBlock: If, localLookup: Map[String, GraphNodePort.OutputPort]): ErrorOr[ConditionalNodeWithNewNodes] = {
    val ifConditionExpression = WdlWomExpression(ifBlock.condition, Option(ifBlock))
    val ifConditionGraphInputExpressionValidation = WdlWomExpression.toExpressionNode(WomIdentifier("conditional"), ifConditionExpression, localLookup, Map.empty)
    val ifConditionTypeValidation = ifConditionExpression.evaluateType(localLookup.map { case (k, v) => k -> v.womType }) flatMap {
      case WdlBooleanType => Valid(())
      case other => s"An if block must be given a boolean expression but instead got '${ifBlock.condition.toWdlString}' (a ${other.toWdlString})".invalidNel
    }

    val innerGraphValidation: ErrorOr[Graph] = WdlGraphNode.buildWomGraph(ifBlock, Set.empty, localLookup)

    (ifConditionGraphInputExpressionValidation, ifConditionTypeValidation, innerGraphValidation) mapN { (ifConditionGraphInputExpression, _, innerGraph) =>
      ConditionalNode.wireInConditional(innerGraph, ifConditionGraphInputExpression)
    }
  }
}
