package wdl.draft3.transforms.ast2wdlom

// TODO 2.11: Either.map doesn't exist in 2.11. When we're 2.12 only, remove this import of cats.syntax.either._
import cats.syntax.either._

import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.model.draft3.elements.{InputDeclarationElement, InputsSectionElement}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._

object AstToInputsSectionElement {

  def convert(a: Ast): Checked[InputsSectionElement] = {

    a.getAttributeAsVector[InputDeclarationElement]("inputs") map { declarations =>
      InputsSectionElement(declarations)
    }

  }
}
