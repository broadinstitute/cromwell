package cromwell.binding

import java.io.File

import cromwell.binding.types._
import cromwell.parser.WdlParser.{Ast, AstNode, Terminal}

case class WdlExpression(ast: AstNode) {
    def evaluate(lookup: String => WdlValue, functions: WdlFunctions): WdlValue = ???
}

trait WdlFunctions {
    type WdlFunction = Seq[WdlExpression] => WdlValue
    def getFunction(name: String): WdlFunction
}
