package wdl4s.wdl

import wdl4s.wdl.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.Ast
import wdl4s.wdl.types.WdlType
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.expression.PlaceholderWomExpression

object TaskOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): TaskOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    TaskOutput(name, wdlType, expression, ast, parent)
  }
  
  def buildWomOutputDefinition(taskOutput: TaskOutput) = OutputDefinition(taskOutput.unqualifiedName, taskOutput.wdlType, PlaceholderWomExpression(Set.empty, taskOutput.wdlType))
}

final case class TaskOutput(unqualifiedName: String, wdlType: WdlType, requiredExpression: WdlExpression, ast: Ast, override val parent: Option[Scope]) extends Output {
  lazy val womOutputDefinition = TaskOutput.buildWomOutputDefinition(this)
}
