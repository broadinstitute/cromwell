package wdl.draft3.transforms.ast2wdlom

// TODO 2.11: Remove import for cats.syntax.either._
import cats.syntax.either._
import common.Checked
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{MetaKvPair, ParameterMetaSectionElement, MetaValueElement}

object AstToParameterMetaSectionElement {
  def convert(ast: Ast): Checked[ParameterMetaSectionElement] =  {
    ast.getAttributeAsVector[MetaKvPair]("map") map { attributes =>
      val asMap = attributes.map(kv => kv.key -> kv.value).toMap
      ParameterMetaSectionElement(asMap)
    }
  }
}
