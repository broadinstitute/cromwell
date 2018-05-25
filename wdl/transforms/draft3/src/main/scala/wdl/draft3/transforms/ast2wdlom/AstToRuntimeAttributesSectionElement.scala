package wdl.draft3.transforms.ast2wdlom

import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements.RuntimeAttributesSectionElement

object AstToRuntimeAttributesSectionElement {
  def convert(ast: Ast): Checked[RuntimeAttributesSectionElement] =  {
    ast.getAttributeAsVector[KvPair]("map") map { attributes =>
      RuntimeAttributesSectionElement(attributes)
    }
  }
}
