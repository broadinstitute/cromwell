package wdl.draft3.transforms.ast2wdlom

import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{OutputDeclarationElement, OutputsSectionElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import common.Checked


object AstToOutputsSectionElement {
  def convert(a: WdlParser.Ast): Checked[OutputsSectionElement] = {
      a.getAttributeAsVector[OutputDeclarationElement]("outputs") map OutputsSectionElement
  }
}
