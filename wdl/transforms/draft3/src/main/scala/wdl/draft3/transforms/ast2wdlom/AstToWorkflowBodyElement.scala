package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.WorkflowBodyElement
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToWorkflowBodyElement {
  def convert(ast: Ast): Checked[WorkflowBodyElement] = ast.getName match {
    case "Inputs" => astNodeToInputsSectionElement(ast)
    case "Outputs" => astNodeToOutputsSectionElement(ast)
    case other => s"No conversion defined for Ast with name $other to WorkflowBodyElement".invalidNelCheck
  }

}
