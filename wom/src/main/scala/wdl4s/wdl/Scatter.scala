package wdl4s.wdl

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
import wdl4s.parser.WdlParser.{Ast, Terminal}
import wdl4s.wdl.types.WdlArrayType
import wdl4s.wom.graph.{Graph, GraphNodePort, RequiredGraphInputNode, ScatterNode}
import wdl4s.wom.graph.ScatterNode.ScatterNodeWithInputs

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

  def womScatterNode(scatter: Scatter, localLookup: Map[String, GraphNodePort.OutputPort]): ErrorOr[ScatterNodeWithInputs] = {
    val scatterCollectionExpression = WdlWomExpression(scatter.collection, Option(scatter))
    val scatterVariableSourceValidation = WdlWomExpression.findInputsforExpression(scatter.item, scatterCollectionExpression, localLookup, Map.empty)
    val scatterItemTypeValidation = scatterCollectionExpression.evaluateType(localLookup.map { case (k, v) => k -> v.womType }) flatMap {
      case WdlArrayType(itemType) => Valid(itemType) // Covers maps because this is a custom unapply
      case other => s"Cannot scatter over a non-array type ${other.toWdlString}".invalidNel
    }

    val innerGraphAndScatterItemInputValidation: ErrorOr[(Graph, RequiredGraphInputNode)] = for {
      scatterItemType <- scatterItemTypeValidation
      womInnerGraphScatterVariableInput = RequiredGraphInputNode(scatter.item, scatterItemType)
      g <- WdlGraphNode.buildWomGraph(scatter, Set(womInnerGraphScatterVariableInput), localLookup)
    } yield (g, womInnerGraphScatterVariableInput)

    (scatterVariableSourceValidation, innerGraphAndScatterItemInputValidation).tupled flatMap { case (scatterVariableSource, (innerGraph, scatterItemInputNode)) =>
      ScatterNode.scatterOverGraph(innerGraph, scatterVariableSource, scatterItemInputNode, localLookup)
    }
  }
}
