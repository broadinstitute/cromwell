package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.ParameterMetaSectionElement
import wom.callable.MetaKvPair

object AstToParameterMetaSectionElement {
  def astToParameterMetaSectionElement(implicit
    astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, MetaKvPair]
  ): CheckedAtoB[GenericAst, ParameterMetaSectionElement] =
    CheckedAtoB.fromCheck("convert AST to parameter_meta section") { ast =>
      ast.getAttributeAsVector[MetaKvPair]("map") map { attributes =>
        val asMap = attributes.map(kv => kv.key -> kv.value).toMap
        ParameterMetaSectionElement(asMap)
      }
    }
}
