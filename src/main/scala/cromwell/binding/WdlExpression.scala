package cromwell.binding

import cromwell.parser.WdlParser.AstNode

case class WdlExpression(ast: AstNode) {
  def evaluate(lookup: String => WdlValue, functions: WdlFunctions): WdlValue = ???
}

trait WdlFunctions {
  type WdlFunction = Seq[WdlExpression] => WdlValue

  def getFunction(name: String): WdlFunction
}
