package wdl.transforms.base.wdlom2wom.expression.renaming

import simulacrum.typeclass
import wdl.model.draft3.elements.ExpressionElement

@typeclass
trait IdentifierLookupRenamer[A <: ExpressionElement] {
  def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): A
}
