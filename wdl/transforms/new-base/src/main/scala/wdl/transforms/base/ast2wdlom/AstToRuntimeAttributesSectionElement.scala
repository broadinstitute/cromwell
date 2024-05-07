package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements.RuntimeAttributesSectionElement

object AstToRuntimeAttributesSectionElement {
  def astToRuntimeSectionElement(implicit
    astNodeToKvPair: CheckedAtoB[GenericAstNode, KvPair]
  ): CheckedAtoB[GenericAst, RuntimeAttributesSectionElement] =
    CheckedAtoB.fromCheck("convert AST to runtime section") { ast =>
      ast.getAttributeAsVector[KvPair]("map") map { attributes =>
        RuntimeAttributesSectionElement(attributes)
      }
    }
}
