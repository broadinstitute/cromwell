package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._

object UnaryOperatorEvaluators {

  implicit val unaryPlusEvaluator: IdentifierLookupRenamer[UnaryPlus] = forOperation(UnaryPlus)
  implicit val unaryNegationEvaluator: IdentifierLookupRenamer[UnaryNegation] = forOperation(UnaryNegation)
  implicit val logicalNotEvaluator: IdentifierLookupRenamer[LogicalNot] = forOperation(LogicalNot)

  private def forOperation[A <: UnaryOperation](constructor: ExpressionElement => A) = new IdentifierLookupRenamer[A] {
    override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
    ): A =
      constructor.apply(expressionElementRenamer.renameIdentifiers(a.argument, renamingMap)(expressionElementRenamer))
  }
}
