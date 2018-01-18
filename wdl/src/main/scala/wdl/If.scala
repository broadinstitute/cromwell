package wdl

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl4s.parser.WdlParser.Ast
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.PlainAnonymousExpressionNode
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

  final lazy val upstreamReferences = condition.variableReferences(this)

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
    val ifConditionExpression = WdlWomExpression(ifBlock.condition, ifBlock)
    val ifConditionGraphInputExpressionValidation = WdlWomExpression.toAnonymousExpressionNode(
      WomIdentifier("conditional"), ifConditionExpression, localLookup, outerLookup, preserveIndexForOuterLookups, ifBlock, PlainAnonymousExpressionNode.apply)
    val ifConditionTypeValidation = ifConditionExpression.evaluateType((localLookup ++ outerLookup).map { case (k, v) => k -> v.womType }) flatMap {
      case WomBooleanType => Valid(())
      case other => s"An if block must be given a boolean expression but instead got '${ifBlock.condition.toWomString}' (a ${other.toDisplayString})".invalidNel
    }

    /**
      * Why? Imagine that we're building three nested levels of a innerGraph.
      * - Say we're building the middle layer.
      * - We have a set of OutputPorts in the outer layer that we can make OGINs to if we need them.
      * - We know that the inner graph might want to make use of those output ports, but we don't know which.
      * - So, we can make OGINs at this layer for all possible OutputPorts in the outer graph and let the inner graph
      * use however many of them it needs.
      */
    val possiblyNeededNestedOgins: Map[String, OuterGraphInputNode] = outerLookup filterNot { case (name, _) => localLookup.contains(name) } map { case (name, outerPort) =>
      name -> OuterGraphInputNode(WomIdentifier(name), outerPort, preserveScatterIndex = preserveIndexForOuterLookups)
    }
    val possiblyNeededNestedOginPorts: Map[String, OutputPort] = possiblyNeededNestedOgins map { case (name: String, ogin: OuterGraphInputNode) => name -> ogin.singleOutputPort }

    val innerGraphValidation: ErrorOr[Graph] = WdlGraphNode.buildWomGraph(
      ifBlock,
      Set.empty,
      outerLookup = localLookup ++ possiblyNeededNestedOginPorts,
      preserveIndexForOuterLookups = true
    )

    (ifConditionGraphInputExpressionValidation, ifConditionTypeValidation, innerGraphValidation) mapN { (ifConditionGraphInputExpression, _, innerGraph) =>
      ConditionalNode.wireInConditional(innerGraph, ifConditionGraphInputExpression)
    }
  }
}
