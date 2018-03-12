package wdl.draft3.transforms.ast2wdlom

import cats.syntax.either._
import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{CommandPartElement, CommandSectionElement}

object AstToCommandSectionElement {
  def convert(ast: Ast): Checked[CommandSectionElement] = {
    ast.getAttributeAsVector[CommandPartElement]("parts") map { parts =>
      CommandSectionElement(parts)
    }
  }

}
