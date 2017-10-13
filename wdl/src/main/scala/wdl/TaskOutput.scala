package wdl

import wdl4s.parser.WdlParser.Ast
import wdl.AstTools.EnhancedAstNode
import wom.callable.Callable.OutputDefinition
import wom.graph._
import wom.types.WdlType

object TaskOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): TaskOutput = {
    val wdlType = ast.getAttribute("type").wdlType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    TaskOutput(name, wdlType, expression, ast, parent)
  }

  def buildWomOutputDefinition(taskOutput: TaskOutput) = {
    OutputDefinition(LocalName(taskOutput.unqualifiedName), taskOutput.wdlType, WdlWomExpression(taskOutput.requiredExpression, from = taskOutput.parent))
  }
}

final case class TaskOutput(unqualifiedName: String, wdlType: WdlType, requiredExpression: WdlExpression, ast: Ast, override val parent: Option[Scope]) extends Output {
  lazy val womOutputDefinition = TaskOutput.buildWomOutputDefinition(this)
}
