package wdl

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl4s.parser.WdlParser.{Ast, Terminal}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.ScatterNode.ScatterNodeWithNewNodes
import wom.graph._
import wom.graph.expression.PlainAnonymousExpressionNode
import wom.types.WomArrayType

/**
  * Scatter class.
  * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
  * @param item Item which this block is scattering over
  * @param collection Wdl Expression corresponding to the collection this scatter is looping through
  */
case class Scatter(index: Int, item: String, collection: WdlExpression, ast: Ast) extends WdlGraphNodeWithUpstreamReferences with WorkflowScoped {
  val unqualifiedName = s"${Scatter.FQNIdentifier}_$index"
  override def appearsInFqn = false

  final lazy val upstreamReferences = collection.variableReferences(this)

  override def toString: String = s"[Scatter fqn=$fullyQualifiedName, item=$item, collection=${collection.toWomString}]"
}

object Scatter {
  val FQNIdentifier = "$scatter"

  /**
    * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
    */
  def apply(ast: Ast, index: Int): Scatter = {
    val item = ast.getAttribute("item").asInstanceOf[Terminal].getSourceString
    new Scatter(index, item, WdlExpression(ast.getAttribute("collection")), ast)
  }

  /**
    * @param preserveIndexForOuterLookups When we're evaluating the scatter collection, should we preserve scatter index when we have to use the outerLookup?
    */
  def womScatterNode(scatter: Scatter, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[ScatterNodeWithNewNodes] = {
    // Convert the scatter collection WdlExpression to a WdlWomExpression 
    val scatterCollectionExpression = WdlWomExpression(scatter.collection, scatter)
    // Generate an ExpressionNode from the WdlWomExpression
    val scatterCollectionExpressionNode =
      WdlWomExpression.toAnonymousExpressionNode(WomIdentifier(scatter.item), scatterCollectionExpression, localLookup, outerLookup, preserveIndexForOuterLookups, scatter, PlainAnonymousExpressionNode.apply)
    // Validate the collection evaluates to a traversable type
    val scatterItemTypeValidation = scatterCollectionExpression.evaluateType((localLookup ++ outerLookup).map { case (k, v) => k -> v.womType }) flatMap {
      case WomArrayType(itemType) => Valid(itemType) // Covers maps because this is a custom unapply (see WdlArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.toDisplayString}".invalidNel
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

    for {
      itemType <- scatterItemTypeValidation
      expressionNode <- scatterCollectionExpressionNode
      // Graph input node for the scatter variable in the inner graph. Note that the type is the array's member type
      womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatter.item), expressionNode.singleExpressionOutputPort, itemType)
      g <- WdlGraphNode.buildWomGraph(scatter, Set(womInnerGraphScatterVariableInput), localLookup ++ possiblyNeededNestedOginPorts, preserveIndexForOuterLookups = false)
    } yield ScatterNode.scatterOverGraph(g, expressionNode, womInnerGraphScatterVariableInput)
  }
}
