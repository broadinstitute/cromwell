package wdl

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl4s.parser.WdlParser.Ast
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.WomBooleanType
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

  override def toString: String = s"[If fqn=$fullyQualifiedName, condition=${condition.toWomString}]"
}

object If {
  val FQNIdentifier = "$if"

  /**
    * @param index Index of the if block. The index is computed during tree generation to reflect WDL scope structure.
    */
  def apply(ast: Ast, index: Int): If = {
    new If(index, WdlExpression(ast.getAttribute("expression")), ast)
  }

  /**
    * @param preserveIndexForOuterLookups When we're evaluating the condition boolean, should we preserve scatter index if we have to use the outerLookup?
    */
  def womConditionalNode(ifBlock: If, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[ConditionalNodeWithNewNodes] = {
    val ifConditionExpression = WdlWomExpression(ifBlock.condition, Option(ifBlock))
    val ifConditionGraphInputExpressionValidation = WdlWomExpression.toExpressionNode(WomIdentifier("conditional"), ifConditionExpression, localLookup, outerLookup, preserveIndexForOuterLookups)
    val ifConditionTypeValidation = ifConditionExpression.evaluateType(localLookup.map { case (k, v) => k -> v.womType }) flatMap {
      case WomBooleanType => Valid(())
      case other => s"An if block must be given a boolean expression but instead got '${ifBlock.condition.toWomString}' (a ${other.toDisplayString})".invalidNel
    }

    val innerGraphValidation: ErrorOr[Graph] = WdlGraphNode.buildWomGraph(
      ifBlock,
      Set.empty,
      outerLookup = localLookup ++ outerLookup,
      preserveIndexForOuterLookups = true
    )

    (ifConditionGraphInputExpressionValidation, ifConditionTypeValidation, innerGraphValidation) mapN { (ifConditionGraphInputExpression, _, innerGraph) =>
      ConditionalNode.wireInConditional(innerGraph, ifConditionGraphInputExpression)
    }
  }
}
