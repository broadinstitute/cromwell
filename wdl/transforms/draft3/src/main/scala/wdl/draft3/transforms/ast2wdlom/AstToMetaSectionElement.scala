package wdl.draft3.transforms.ast2wdlom

import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.MetaSectionElement
import wom.callable.MetaKvPair

object AstToMetaSectionElement {
  def convert(ast: Ast): Checked[MetaSectionElement] =  {
    ast.getAttributeAsVector[MetaKvPair]("map") map { attributes =>
      val asMap = attributes.map(kv => kv.key -> kv.value).toMap
      MetaSectionElement(asMap)
    }
  }
}
