package wdl

import cats.data.Validated.Valid
import cats.syntax.validated._
import lenthall.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.types.WdlArrayType
import wdl4s.parser.WdlParser.{Ast, Terminal}
import wom.graph.ScatterNode.ScatterNodeWithNewNodes
import wom.graph._

/**
  * Scatter class.
  * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
  * @param item Item which this block is scattering over
  * @param collection Wdl Expression corresponding to the collection this scatter is looping through
  */
case class Scatter(index: Int, item: String, collection: WdlExpression, ast: Ast) extends WdlGraphNodeWithUpstreamReferences with WorkflowScoped {
  val unqualifiedName = s"${Scatter.FQNIdentifier}_$index"
  override def appearsInFqn = false

  final val upstreamReferences = collection.variableReferences

  override def toString: String = s"[Scatter fqn=$fullyQualifiedName, item=$item, collection=${collection.toWdlString}]"
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

  def womScatterNode(scatter: Scatter, localLookup: Map[String, GraphNodePort.OutputPort]): ErrorOr[ScatterNodeWithNewNodes] = {
    // Convert the scatter collection WdlExpression to a WdlWomExpression 
    val scatterCollectionExpression = WdlWomExpression(scatter.collection, Option(scatter))
    // Generate an ExpressionNode from the WdlWomExpression
    val scatterCollectionExpressionNode = WdlWomExpression.toExpressionNode(WomIdentifier(scatter.item), scatterCollectionExpression, localLookup, Map.empty)
    // Validate the collection evaluates to a traversable type
    val scatterCollectionTypeValidation = scatterCollectionExpression.evaluateType(localLookup.map { case (k, v) => k -> v.womType }) flatMap {
      case collectionType: WdlArrayType => Valid(collectionType) // Covers maps because this is a custom unapply (see WdlArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.toWdlString}".invalidNel
    }

    for {
      collectionType <- scatterCollectionTypeValidation
      expressionNode <- scatterCollectionExpressionNode
      // Graph input node for the scatter variable in the inner graph. Note that the type is the array's member type
      womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatter.item), expressionNode.singleExpressionOutputPort, collectionType.memberType)
      g <- WdlGraphNode.buildWomGraph(scatter, Set(womInnerGraphScatterVariableInput), localLookup)
    } yield ScatterNode.scatterOverGraph(g, expressionNode, womInnerGraphScatterVariableInput)
  }
}
