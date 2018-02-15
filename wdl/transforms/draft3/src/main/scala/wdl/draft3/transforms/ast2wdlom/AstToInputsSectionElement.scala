package wdl.draft3.transforms.ast2wdlom

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
