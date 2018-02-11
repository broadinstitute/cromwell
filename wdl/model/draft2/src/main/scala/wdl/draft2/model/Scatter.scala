package wdl.draft2.model

import wdl.draft2.parser.WdlParser.{Ast, Terminal}

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
}
