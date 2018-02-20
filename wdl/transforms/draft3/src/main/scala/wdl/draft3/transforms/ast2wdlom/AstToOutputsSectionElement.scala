package wdl.draft3.transforms.ast2wdlom

// TODO 2.11: Either.map doesn't exist in 2.11. When we're 2.12 only, remove this import of cats.syntax.either._
import cats.syntax.either._

import wdl.draft3.parser.WdlParser
import wdl.model.draft3.elements.{DeclarationElement, OutputsSectionElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import common.Checked

object AstToOutputsSectionElement {
  def convert(a: WdlParser.Ast): Checked[OutputsSectionElement] = {
      a.getAttributeAsVector[DeclarationElement]("outputs") map OutputsSectionElement
  }
}
