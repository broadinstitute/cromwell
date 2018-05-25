package wdl.draft3.transforms.ast2wdlom

import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{DeclarationContent, OutputDeclarationElement, OutputsSectionElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import common.Checked

object AstToOutputsSectionElement {
  def convert(a: WdlParser.Ast): Checked[OutputsSectionElement] =
      a.getAttributeAsVector[DeclarationContent]("outputs").map(_.map(OutputDeclarationElement.fromContent)).map(OutputsSectionElement)
}
