package wdl.draft3.transforms.ast2wdlom

// TODO 2.11: Remove this import to cats.syntax.either._
import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{IntermediateValueDeclarationElement, WorkflowGraphElement}

object AstToWorkflowGraphNodeElement {
  def convert(ast: Ast): Checked[WorkflowGraphElement] = ast.getName match {
    case "Declaration" => astNodeToDeclarationContent(ast).map(IntermediateValueDeclarationElement.fromContent)
    case "Call" => astNodeToCallElement(ast)
    case "Scatter" => astNodeToScatterElement(ast)
    case "If" => astNodeToIfElement(ast)
    case other => s"No conversion defined for Ast with name $other to WorkflowGraphElement".invalidNelCheck
  }
}
