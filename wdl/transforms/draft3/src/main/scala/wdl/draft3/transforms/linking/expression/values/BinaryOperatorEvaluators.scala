package wdl.draft3.transforms.linking.expression.values

import common.validation.Validation._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.IoFunctionSet
import wom.values.WomValue

import scala.util.Try


object BinaryOperatorEvaluators {
  implicit val logicalOrEvaluator: ValueEvaluator[LogicalOr] = forOperation(_.or(_))
  implicit val logicalAndEvaluator: ValueEvaluator[LogicalAnd] = forOperation(_.and(_))
  implicit val equalsEvaluator: ValueEvaluator[Equals] = forOperation(_.equals(_))
  implicit val notEqualsEvaluator: ValueEvaluator[NotEquals] = forOperation(_.notEquals(_))
  implicit val lessThanEvaluator: ValueEvaluator[LessThan] = forOperation(_.lessThan(_))
  implicit val lessThanOrEqualEvaluator: ValueEvaluator[LessThanOrEquals] = forOperation(_.lessThanOrEqual(_))
  implicit val greaterThanEvaluator: ValueEvaluator[GreaterThan] = forOperation(_.greaterThan(_))
  implicit val greaterThanOrEqualEvaluator: ValueEvaluator[GreaterThanOrEquals] = forOperation(_.greaterThanOrEqual(_))
  implicit val addEvaluator: ValueEvaluator[Add] = forOperation(_.add(_))
  implicit val subtractEvaluator: ValueEvaluator[Subtract] = forOperation(_.subtract(_))
  implicit val multiplyEvaluator: ValueEvaluator[Multiply] = forOperation(_.multiply(_))
  implicit val divideEvaluator: ValueEvaluator[Divide] = forOperation(_.divide(_))
  implicit val remainderEvaluator: ValueEvaluator[Remainder] = forOperation(_.mod(_))

  private def forOperation[A <: BinaryOperation](op: (WomValue, WomValue) => Try[WomValue]) = new ValueEvaluator[A] {
    override def evaluateValue(a: A,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      (a.left.evaluateValue(inputs, ioFunctionSet),
        a.right.evaluateValue(inputs, ioFunctionSet)) flatMapN { (left, right) => op(left, right).toErrorOr }
  }
}
