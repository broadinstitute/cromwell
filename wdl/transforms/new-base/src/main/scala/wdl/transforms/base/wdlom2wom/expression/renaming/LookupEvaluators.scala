package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._

object LookupEvaluators {

  implicit val identifierLookupEvaluator: IdentifierLookupRenamer[IdentifierLookup] =
    new IdentifierLookupRenamer[IdentifierLookup] {
      override def renameIdentifiers(a: IdentifierLookup, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): IdentifierLookup = IdentifierLookup(
        if (renamingMap.contains(a.identifier))
          renamingMap(a.identifier)
        else
          a.identifier
      )
    }

  implicit val expressionMemberAccessEvaluator: IdentifierLookupRenamer[ExpressionMemberAccess] =
    new IdentifierLookupRenamer[ExpressionMemberAccess] {
      override def renameIdentifiers(a: ExpressionMemberAccess, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): ExpressionMemberAccess = ExpressionMemberAccess(
        expressionElementRenamer.renameIdentifiers(a.expression, renamingMap)(expressionElementRenamer),
        a.memberAccessTail
      )
    }

  implicit val identifierMemberAccessEvaluator: IdentifierLookupRenamer[IdentifierMemberAccess] =
    new IdentifierLookupRenamer[IdentifierMemberAccess] {
      override def renameIdentifiers(a: IdentifierMemberAccess, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): IdentifierMemberAccess = IdentifierMemberAccess(
        if (renamingMap.contains(a.first)) renamingMap(a.first) else a.first,
        a.second,
        a.memberAccessTail
      )
    }

  implicit val indexAccessIdentifierLookupRenamer: IdentifierLookupRenamer[IndexAccess] =
    new IdentifierLookupRenamer[IndexAccess] {
      override def renameIdentifiers(a: IndexAccess, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): IndexAccess = IndexAccess(
        expressionElementRenamer.renameIdentifiers(a.expressionElement, renamingMap)(expressionElementRenamer),
        expressionElementRenamer.renameIdentifiers(a.index, renamingMap)(expressionElementRenamer)
      )
    }
}
