package wdl.draft2.model

import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.parser.WdlParser.Ast
import wom.types.WomType

object WorkflowOutput {
  def apply(ast: Ast, syntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): WorkflowOutput = {
    val womType = ast.getAttribute("type").womType(syntaxErrorFormatter)
    val name = ast.getAttribute("name").sourceString
    val expression = WdlExpression(ast.getAttribute("expression"))
    WorkflowOutput(name, womType, expression, ast, parent)
  }
}

case class WorkflowOutput(unqualifiedName: String,
                          womType: WomType,
                          requiredExpression: WdlExpression,
                          ast: Ast,
                          override val parent: Option[Scope]
) extends Output
