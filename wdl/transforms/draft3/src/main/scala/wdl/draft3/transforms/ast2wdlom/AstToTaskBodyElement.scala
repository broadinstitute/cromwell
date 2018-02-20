package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{IntermediateValueDeclarationElement, TaskBodyElement}

object AstToTaskBodyElement {
  def convert(ast: Ast): Checked[TaskBodyElement] = ast.getName match {
    case "Inputs" => astNodeToInputsSectionElement(ast)
    case "Outputs" => astNodeToOutputsSectionElement(ast)
    case "Declaration" => astNodeToDeclarationContent(ast).map(IntermediateValueDeclarationElement.fromContent)
    case other => s"No conversion defined for Ast with name $other to TaskBodyElement".invalidNelCheck
  }
}
