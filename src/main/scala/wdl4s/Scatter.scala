package wdl4s

import wdl4s.parser.WdlParser.{Ast, Terminal}
import wdl4s.AstTools.EnhancedAstNode

object Scatter {
  val FQNIdentifier = "$scatter"

  /**
   * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
   */
  def apply(ast: Ast, index: Int): Scatter = {
    val item = ast.getAttribute("item").asInstanceOf[Terminal].getSourceString
    new Scatter(index, item, WdlExpression(ast.getAttribute("collection")), ast)
  }
}

/**
 * Scatter class.
 * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
 * @param item Item which this block is scattering over
 * @param collection Wdl Expression corresponding to the collection this scatter is looping through
 */
case class Scatter(index: Int, item: String, collection: WdlExpression, ast: Ast) extends GraphNodeWithUpstreamReferences with WorkflowScoped {
  val unqualifiedName = s"${Scatter.FQNIdentifier}_$index"
  override def appearsInFqn = false

  final val upstreamReferences = collection.variableReferences

  override def toString: String = s"[Scatter fqn=$fullyQualifiedName, item=$item, collection=${collection.toWdlString}]"
}
