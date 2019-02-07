package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._


object LiteralEvaluators {
  implicit val primitiveIdentifierLookupRenamer: IdentifierLookupRenamer[PrimitiveLiteralExpressionElement] = new IdentifierLookupRenamer[PrimitiveLiteralExpressionElement] {
    override def renameIdentifiers(a: PrimitiveLiteralExpressionElement, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): PrimitiveLiteralExpressionElement = a
  }

  implicit val stringLiteralEvaluator: IdentifierLookupRenamer[StringLiteral] = new IdentifierLookupRenamer[StringLiteral] {
    override def renameIdentifiers(a: StringLiteral, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): StringLiteral = a
  }

  implicit val stringExpressionEvaluator: IdentifierLookupRenamer[StringExpression] = new IdentifierLookupRenamer[StringExpression] {
    override def renameIdentifiers(a: StringExpression, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): StringExpression = StringExpression(
      a.pieces map {
        case sp: StringPlaceholder => StringPlaceholder(expressionElementRenamer.renameIdentifiers(sp.expr, renamingMap)(expressionElementRenamer))
        case other => other
      }
    )
  }

  implicit val objectLiteralEvaluator: IdentifierLookupRenamer[ObjectLiteral] = new IdentifierLookupRenamer[ObjectLiteral] {
    override def renameIdentifiers(a: ObjectLiteral, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): ObjectLiteral = ObjectLiteral(
      a.elements map {
        case (key, value) => key -> expressionElementRenamer.renameIdentifiers(value, renamingMap)(expressionElementRenamer)
      }
    )
  }

  implicit val mapLiteralEvaluator: IdentifierLookupRenamer[MapLiteral] = new IdentifierLookupRenamer[MapLiteral] {
    override def renameIdentifiers(a: MapLiteral, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): MapLiteral = MapLiteral(
      a.elements map {
        case (key, value) => expressionElementRenamer.renameIdentifiers(key, renamingMap)(expressionElementRenamer) -> expressionElementRenamer.renameIdentifiers(value, renamingMap)(expressionElementRenamer)
      }
    )
  }

  implicit val arrayLiteralEvaluator: IdentifierLookupRenamer[ArrayLiteral] = new IdentifierLookupRenamer[ArrayLiteral] {
    override def renameIdentifiers(a: ArrayLiteral, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): ArrayLiteral = ArrayLiteral(
      a.elements map { e => expressionElementRenamer.renameIdentifiers(e, renamingMap)(expressionElementRenamer) }
    )
  }

  implicit val pairLiteralEvaluator: IdentifierLookupRenamer[PairLiteral] = new IdentifierLookupRenamer[PairLiteral] {
    override def renameIdentifiers(a: PairLiteral, renamingMap: Map[String, String])
                                  (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]): PairLiteral = PairLiteral(
      expressionElementRenamer.renameIdentifiers(a.left, renamingMap)(expressionElementRenamer),
      expressionElementRenamer.renameIdentifiers(a.right, renamingMap)(expressionElementRenamer)
    )
  }
}
