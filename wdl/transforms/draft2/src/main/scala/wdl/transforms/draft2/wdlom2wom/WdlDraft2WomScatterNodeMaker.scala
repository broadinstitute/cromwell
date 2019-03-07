package wdl.transforms.draft2.wdlom2wom

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.draft2.model.{Scatter, Scope, WdlWomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.ScatterNode.ScatterNodeWithNewNodes
import wom.graph._
import wom.graph.expression.PlainAnonymousExpressionNode
import wom.transforms.WomScatterNodeMaker
import wom.transforms.WomGraphMaker.ops._
import wom.types.WomArrayType

object WdlDraft2WomScatterNodeMaker extends WomScatterNodeMaker[Scatter] {
  override def toWomScatterNode(scatter: Scatter,
                                localLookup: Map[String, GraphNodePort.OutputPort],
                                outerLookup: Map[String, GraphNodePort.OutputPort],
                                preserveIndexForOuterLookups: Boolean,
                                inASubworkflow: Boolean): ErrorOr[ScatterNodeWithNewNodes] = {

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
        * preserveScatterIndex = false because in the absence of support of nested scatters,
        * the index should never be preserved when for nodes coming from outside the scatter.
       */
      name -> OuterGraphInputNode(WomIdentifier(name), outerPort, preserveScatterIndex = false)
    }
    val possiblyNeededNestedOginPorts: Map[String, OutputPort] = possiblyNeededNestedOgins map { case (name: String, ogin: OuterGraphInputNode) => name -> ogin.singleOutputPort }

    // Convert the scatter collection WdlExpression to a WdlWomExpression
    val scatterCollectionExpression = WdlWomExpression(scatter.collection, scatter)
    // Generate an ExpressionNode from the WdlWomExpression
    val scatterCollectionExpressionNode =
      WdlWomExpression.toAnonymousExpressionNode(WomIdentifier(scatter.item), scatterCollectionExpression, localLookup ++ possiblyNeededNestedOginPorts, Map.empty, preserveIndexForOuterLookups, scatter, PlainAnonymousExpressionNode.apply)
    // Validate the collection evaluates to a traversable type
    val scatterItemTypeValidation = scatterCollectionExpression.evaluateType((localLookup ++ outerLookup).map { case (k, v) => k -> v.womType }) flatMap {
      case WomArrayType(itemType) => Valid(itemType) // Covers maps because this is a custom unapply (see WdlArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.stableName}".invalidNel
    }

    for {
      itemType <- scatterItemTypeValidation
      expressionNode <- scatterCollectionExpressionNode
      // Graph input node for the scatter variable in the inner graph. Note that the type is the array's member type
      womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatter.item), expressionNode, itemType)
      g <- (scatter: Scope).toWomGraph(Set(womInnerGraphScatterVariableInput), localLookup ++ possiblyNeededNestedOginPorts, preserveIndexForOuterLookups = false, inASubworkflow)
    } yield ScatterNode.scatterOverGraph(g, womInnerGraphScatterVariableInput)

  }
}
