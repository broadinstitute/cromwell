package wdl.draft3.transforms.ast2wdlom

// TODO 2.11: Remove this import to cats.syntax.either._
import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{IntermediateValueDeclarationElement, WorkflowBodyElement}

object AstToWorkflowBodyElement {
  def convert(ast: Ast): Checked[WorkflowBodyElement] = ast.getName match {
    case "Inputs" => astNodeToInputsSectionElement(ast)
    case "Outputs" => astNodeToOutputsSectionElement(ast)
    case "Declaration" => astNodeToDeclarationContent(ast).map(IntermediateValueDeclarationElement.fromContent)
    case "Call" => astNodeToCallElement(ast)
    case "Meta" => astNodeToMetaSectionElement(ast)
    case "ParameterMeta" => astNodeToParameterMetaSectionElement(ast)
    case other => s"No conversion defined for Ast with name $other to WorkflowBodyElement".invalidNelCheck
  }
}
