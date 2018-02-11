package wdl.draft2.model

import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.parser.WdlParser.Ast
import wom.callable.Callable.OutputDefinition
import wom.graph._
import wom.types.WomType

object TaskOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): TaskOutput = {
    val womType = ast.getAttribute("type").womType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    TaskOutput(name, womType, expression, ast, parent)
  }

  def buildWomOutputDefinition(taskOutput: TaskOutput) = {
    OutputDefinition(LocalName(taskOutput.unqualifiedName), taskOutput.womType, WdlWomExpression(taskOutput.requiredExpression, from = taskOutput))
  }
}

final case class TaskOutput(unqualifiedName: String, womType: WomType, requiredExpression: WdlExpression, ast: Ast, override val parent: Option[Scope]) extends Output {
  lazy val womOutputDefinition: OutputDefinition = TaskOutput.buildWomOutputDefinition(this)
}
