package cromwell.binding

import cromwell.parser.WdlParser.{Ast, Terminal}

object Scatter {
  /**
   * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
   */
  def apply(ast: Ast, index: Int, parent: Option[Scope]): Scatter = {
    val item = ast.getAttribute("item").asInstanceOf[Terminal].getSourceString
    new Scatter(index, item, WdlExpression(ast.getAttribute("collection")), parent)
  }
}

/**
 * Scatter class.
 * @param index Index of the scatter block. The index is computed during tree generation to reflect wdl scatter blocks structure.
 * @param item Item which this block is scattering over
 * @param collection Wdl Expression corresponding to the collection this scatter is looping through
 * @param parent Parent of this scatter
 */
case class Scatter(index: Int, item: String, collection: WdlExpression, parent: Option[Scope]) extends Scope {

  val name = s"$$scatter_$index"

  override def appearsInFQN = false

}
