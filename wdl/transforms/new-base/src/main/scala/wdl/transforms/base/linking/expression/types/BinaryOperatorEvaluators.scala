package wdl.transforms.base.linking.expression.types

import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import common.validation.Validation._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types.WomType

import scala.util.Try

object BinaryOperatorEvaluators {
  implicit val logicalOrEvaluator: TypeEvaluator[LogicalOr] = forOperation(_.or(_))
  implicit val logicalAndEvaluator: TypeEvaluator[LogicalAnd] = forOperation(_.and(_))
  implicit val equalsEvaluator: TypeEvaluator[Equals] = forOperation(_.equalsType(_))
  implicit val notEqualsEvaluator: TypeEvaluator[NotEquals] = forOperation(_.notEquals(_))
  implicit val lessThanEvaluator: TypeEvaluator[LessThan] = forOperation(_.lessThan(_))
  implicit val lessThanOrEqualEvaluator: TypeEvaluator[LessThanOrEquals] = forOperation(_.lessThanOrEqual(_))
  implicit val greaterThanEvaluator: TypeEvaluator[GreaterThan] = forOperation(_.greaterThan(_))
  implicit val greaterThanOrEqualEvaluator: TypeEvaluator[GreaterThanOrEquals] = forOperation(_.greaterThanOrEqual(_))
  implicit val addEvaluator: TypeEvaluator[Add] = forOperation(_.add(_))
  implicit val subtractEvaluator: TypeEvaluator[Subtract] = forOperation(_.subtract(_))
  implicit val multiplyEvaluator: TypeEvaluator[Multiply] = forOperation(_.multiply(_))
  implicit val divideEvaluator: TypeEvaluator[Divide] = forOperation(_.divide(_))
  implicit val remainderEvaluator: TypeEvaluator[Remainder] = forOperation(_.mod(_))

  private def forOperation[A <: BinaryOperation](op: (WomType, WomType) => Try[WomType]) = new TypeEvaluator[A] {
    override def evaluateType(a: A,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      (a.left.evaluateType(linkedValues, typeAliases), a.right.evaluateType(linkedValues, typeAliases)) flatMapN {
        (left, right) =>
          op(left, right).toErrorOr
      }
  }
}
