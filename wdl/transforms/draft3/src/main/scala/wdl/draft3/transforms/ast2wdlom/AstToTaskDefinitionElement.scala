package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.TaskDefinitionElement

object AstToTaskDefinitionElement {

  def convert(a: Ast): Checked[TaskDefinitionElement] =
    "FromAst[TaskDefinitionElement](a: Ast) is not implemented".invalidNelCheck
}
