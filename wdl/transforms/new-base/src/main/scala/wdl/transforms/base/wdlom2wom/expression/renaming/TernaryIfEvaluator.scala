package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._

object TernaryIfEvaluator {
  implicit val ternaryIfEvaluator: IdentifierLookupRenamer[TernaryIf] = new IdentifierLookupRenamer[TernaryIf] {
    override def renameIdentifiers(a: TernaryIf, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): TernaryIf = TernaryIf(
      expressionElementRenamer.renameIdentifiers(a.condition, renamingMap)(expressionElementRenamer),
      expressionElementRenamer.renameIdentifiers(a.ifTrue, renamingMap)(expressionElementRenamer),
      expressionElementRenamer.renameIdentifiers(a.ifFalse, renamingMap)(expressionElementRenamer)
    )
  }
}
