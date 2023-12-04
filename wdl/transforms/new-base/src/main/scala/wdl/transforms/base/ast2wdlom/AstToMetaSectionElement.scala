package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.MetaSectionElement
import wom.callable.MetaKvPair

object AstToMetaSectionElement {
  def astToMetaSectionElement(implicit
    astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, MetaKvPair]
  ): CheckedAtoB[GenericAst, MetaSectionElement] = CheckedAtoB.fromCheck("convert AST to Meta Section") { ast =>
    ast.getAttributeAsVector[MetaKvPair]("map") map { attributes =>
      val asMap = attributes.map(kv => kv.key -> kv.value).toMap
      MetaSectionElement(asMap)
    }
  }
}
