package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.{DeclarationContent, OutputDeclarationElement, OutputsSectionElement}

object AstToOutputsSectionElement {
  def astToOutputSectionElement(implicit astNodeToMetaKvPair: CheckedAtoB[GenericAstNode, DeclarationContent]
                             ): CheckedAtoB[GenericAst, OutputsSectionElement] = CheckedAtoB.fromCheck("read outputs section") { a =>

    a.getAttributeAsVector[DeclarationContent]("outputs").map(_.map(OutputDeclarationElement.fromContent)).map(OutputsSectionElement)
  }
}
