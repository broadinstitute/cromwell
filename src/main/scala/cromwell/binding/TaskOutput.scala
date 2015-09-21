package cromwell.binding

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.expression.{NoFunctions, WdlStandardLibraryFunctionsType}
import cromwell.binding.types.WdlType
import cromwell.parser.WdlParser.{Ast, SyntaxError, Terminal}
import com.typesafe.scalalogging.LazyLogging

import scala.util.{Failure, Success}

object TaskOutput extends LazyLogging {
  def apply(ast: Ast, declarations: Seq[Declaration], syntaxErrorFormatter: WdlSyntaxErrorFormatter): TaskOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("var").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))

    // .get below because contract with the lookup() function is that it signals
    // an error by throwing an exception which is wrapped in a Try() in the
    // evaluator code
    def lookup(n: String) = declarations.find(_.name == n).get.wdlType

    def badCoercionException(expressionWdlType: WdlType) = throw new SyntaxError(syntaxErrorFormatter.taskOutputExpressionTypeDoesNotMatchDeclaredType(
      ast.getAttribute("var").asInstanceOf[Terminal],
      wdlType,
      expressionWdlType
    ))

    expression.evaluateType(lookup, new WdlStandardLibraryFunctionsType) match {
      case Success(expressionWdlType) if !wdlType.isCoerceableFrom(expressionWdlType) => badCoercionException(expressionWdlType)
      case Success(expressionWdlType) =>
        val expressionValue = expression.evaluate((s: String) => throw new Throwable("not implemented"), new NoFunctions)
        expressionValue match {
          case Success(value) if wdlType.coerceRawValue(value).isFailure => badCoercionException(expressionWdlType)
          case _ =>
        }
      case Failure(ex) => logger.error(s"Could not determine type of expression: ${expression.toWdlString}")
    }
    new TaskOutput(name, wdlType, expression)
  }
}

case class TaskOutput(name: String, wdlType: WdlType, expression: WdlExpression)
