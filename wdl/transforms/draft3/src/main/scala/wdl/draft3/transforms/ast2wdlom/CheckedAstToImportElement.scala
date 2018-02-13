package wdl.draft3.transforms.ast2wdlom

import common.Checked
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.ImportElement

object CheckedAstToImportElement {

  type CheckedAstToImportElement = CheckedAstTo[ImportElement]
  def instance: CheckedAstToImportElement = CheckedAtoB(convert _)

  def convert(a: Ast): Checked[ImportElement] =
    "FromAst[ImportElement](a: Ast) is not implemented".invalidNelCheck
}
