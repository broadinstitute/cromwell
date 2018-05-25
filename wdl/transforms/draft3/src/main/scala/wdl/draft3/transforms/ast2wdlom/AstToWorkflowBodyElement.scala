package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.WorkflowBodyElement

object AstToWorkflowBodyElement {
  def convert(ast: Ast): Checked[WorkflowBodyElement] = ast.getName match {
    case "Inputs" => astNodeToInputsSectionElement(ast)
    case "Outputs" => astNodeToOutputsSectionElement(ast)
    case "Meta" => astNodeToMetaSectionElement(ast)
    case "ParameterMeta" => astNodeToParameterMetaSectionElement(ast)
    case "Declaration" | "Call" | "Scatter" | "If" => astNodeToGraphElement(ast)
    case other => s"No conversion defined for Ast with name $other to WorkflowBodyElement".invalidNelCheck
  }
}
