package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctionsType}
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import com.typesafe.scalalogging.LazyLogging

import scala.util.{Failure, Success}

object TaskOutput extends LazyLogging {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter): TaskOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("var").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    TaskOutput(name, wdlType, expression)
  }
}

case class TaskOutput(name: String, wdlType: WdlType, requiredExpression: WdlExpression) extends Declaration {
  override val postfixQuantifier = None
  override val expression = Option(requiredExpression)
}
