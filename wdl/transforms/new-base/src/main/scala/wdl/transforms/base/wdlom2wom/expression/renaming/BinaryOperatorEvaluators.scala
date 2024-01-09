package wdl.transforms.base.wdlom2wom.expression.renaming

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._

object BinaryOperatorEvaluators {

  implicit val logicalOrEvaluator: IdentifierLookupRenamer[LogicalOr] = forOperation(LogicalOr)
  implicit val logicalAndEvaluator: IdentifierLookupRenamer[LogicalAnd] = forOperation(LogicalAnd)
  implicit val equalsEvaluator: IdentifierLookupRenamer[Equals] = forOperation(Equals)
  implicit val notEqualsEvaluator: IdentifierLookupRenamer[NotEquals] = forOperation(NotEquals)
  implicit val lessThanEvaluator: IdentifierLookupRenamer[LessThan] = forOperation(LessThan)
  implicit val lessThanOrEqualEvaluator: IdentifierLookupRenamer[LessThanOrEquals] = forOperation(LessThanOrEquals)
  implicit val greaterThanEvaluator: IdentifierLookupRenamer[GreaterThan] = forOperation(GreaterThan)
  implicit val greaterThanOrEqualEvaluator: IdentifierLookupRenamer[GreaterThanOrEquals] = forOperation(
    GreaterThanOrEquals
  )
  implicit val addEvaluator: IdentifierLookupRenamer[Add] = forOperation(Add)
  implicit val subtractEvaluator: IdentifierLookupRenamer[Subtract] = forOperation(Subtract)
  implicit val multiplyEvaluator: IdentifierLookupRenamer[Multiply] = forOperation(Multiply)
  implicit val divideEvaluator: IdentifierLookupRenamer[Divide] = forOperation(Divide)
  implicit val remainderEvaluator: IdentifierLookupRenamer[Remainder] = forOperation(Remainder)

  private def forOperation[A <: BinaryOperation](constructor: (ExpressionElement, ExpressionElement) => A) =
    new IdentifierLookupRenamer[A] {
      override def renameIdentifiers(a: A, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): A =
        constructor.apply(
          expressionElementRenamer.renameIdentifiers(a.left, renamingMap)(expressionElementRenamer),
          expressionElementRenamer.renameIdentifiers(a.right, renamingMap)(expressionElementRenamer)
        )
    }
}
