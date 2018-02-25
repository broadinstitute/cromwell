package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.FileBodyElement

object AstToFileBodyElement {
  def convert(ast: Ast): Checked[FileBodyElement] = ast.getName match {
    case "Workflow" => astNodeToWorkflowDefinitionElement(ast)
    case "Task" => astNodeToTaskDefinitionElement(ast)
    case "Struct" => astNodeToStructEntry(ast)
    case other => s"No conversion defined for Ast with name $other".invalidNelCheck
  }
}
