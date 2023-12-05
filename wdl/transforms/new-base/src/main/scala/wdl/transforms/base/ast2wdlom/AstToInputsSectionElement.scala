package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.{InputDeclarationElement, InputsSectionElement}

object AstToInputsSectionElement {

  def astToInputsSectionElement(implicit
    astNodeToInputDeclaration: CheckedAtoB[GenericAstNode, InputDeclarationElement]
  ): CheckedAtoB[GenericAst, InputsSectionElement] = CheckedAtoB.fromCheck("convert Ast to Inputs Section") { a =>
    a.getAttributeAsVector[InputDeclarationElement]("inputs") map { declarations =>
      InputsSectionElement(declarations)
    }

  }
}
