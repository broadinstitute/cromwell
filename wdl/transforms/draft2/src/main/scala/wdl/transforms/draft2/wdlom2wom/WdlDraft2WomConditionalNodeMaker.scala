package wdl.transforms.draft2.wdlom2wom

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.draft2.model.{If, Scope, WdlWomExpression}
import wdl.draft2.model.{Scope, WdlWomExpression}
import wom.graph.ConditionalNode.ConditionalNodeWithNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.PlainAnonymousExpressionNode
import wom.transforms.WomConditionalNodeMaker
import wom.transforms.WomGraphMaker.ops._
import wom.types.WomBooleanType

object WdlDraft2WomConditionalNodeMaker extends WomConditionalNodeMaker[If] {
  /**
    * @param preserveIndexForOuterLookups When we're evaluating the condition boolean, should we preserve scatter index if we have to use the outerLookup?
    */
  override def toWomConditionalNode(ifBlock: If, localLookup: Map[String, OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean, inASubworkflow: Boolean): ErrorOr[ConditionalNodeWithNewNodes] = {

    /*
      * Why? Imagine that we're building three nested levels of a innerGraph.
      * - Say we're building the middle layer.
      * - We have a set of OutputPorts in the outer layer that we can make OGINs to if we need them.
      * - We know that the inner graph might want to make use of those output ports, but we don't know which.
      * - So, we can make OGINs at this layer for all possible OutputPorts in the outer graph and let the inner graph
      * use however many of them it needs.
      */
    val possiblyNeededNestedOgins: Map[String, OuterGraphInputNode] = outerLookup filterNot { case (name, _) => localLookup.contains(name) } map { case (name, outerPort) =>
      /*
        * preserveIndexForOuterLookups indicates us whether or not nodes in the outerLookup are in the same scatter inn graph as this node
        * preserveIndexForOuterLookups = false means the outerLookup nodes are outside a scatter containing this conditional node
        * preserveIndexForOuterLookups = true means the above predicate does not hold
        *
        * When creating OGINs from those outer lookup nodes for the inner graph we want to make sure we set their preserveScatterIndex to preserveIndexForOuterLookups
        * because they effectively represent the outer lookup nodes inside the conditional. So whether the index must be preserved depends on whether this
        * conditional node has been asked to "preserveIndexForOuterLookups".
       */
      name -> OuterGraphInputNode(WomIdentifier(name), outerPort, preserveScatterIndex = preserveIndexForOuterLookups)
    }
    val possiblyNeededNestedOginPorts: Map[String, OutputPort] = possiblyNeededNestedOgins map { case (name: String, ogin: OuterGraphInputNode) => name -> ogin.singleOutputPort }

    val ifConditionExpression = WdlWomExpression(ifBlock.condition, ifBlock)
    val ifConditionGraphInputExpressionValidation = WdlWomExpression.toAnonymousExpressionNode(
      WomIdentifier("conditional"), ifConditionExpression, localLookup ++ possiblyNeededNestedOginPorts, Map.empty, preserveIndexForOuterLookups, ifBlock, PlainAnonymousExpressionNode.apply)
    val ifConditionTypeValidation = ifConditionExpression.evaluateType((localLookup ++ outerLookup).map { case (k, v) => k -> v.womType }) flatMap {
      case coerceable if WomBooleanType.isCoerceableFrom(coerceable) => Valid(())
      case other => s"An if block must be given a boolean expression but instead got '${ifBlock.condition.toWomString}' (a ${other.stableName})".invalidNel
    }

    val innerGraphValidation: ErrorOr[Graph] = (ifBlock: Scope).toWomGraph(
      Set.empty,
      outerLookup = localLookup ++ possiblyNeededNestedOginPorts,
      preserveIndexForOuterLookups = true,
      inASubworkflow = inASubworkflow
    )

    (ifConditionGraphInputExpressionValidation, ifConditionTypeValidation, innerGraphValidation) mapN { (ifConditionGraphInputExpression, _, innerGraph) =>
      ConditionalNode.wireInConditional(innerGraph, ifConditionGraphInputExpression)
    }
  }
}
