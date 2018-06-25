package wdl.draft3.transforms.linking.expression.values

import cats.syntax.validated._
import common.validation.Validation._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.OptionalNotSuppliedException
import wom.expression.IoFunctionSet
import wom.values.{WomBoolean, WomString, WomValue}

import scala.util.Try

object BinaryOperatorEvaluators {
  implicit val logicalOrEvaluator: ValueEvaluator[LogicalOr] = forOperationWithShortCircuit(_.or(_), shortCircuit = { case WomBoolean(true) => WomBoolean(true) })
  implicit val logicalAndEvaluator: ValueEvaluator[LogicalAnd] = forOperationWithShortCircuit(_.and(_), shortCircuit = { case WomBoolean(false) => WomBoolean(false) })
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

  private def forOperation[A <: BinaryOperation](op: (WomValue, WomValue) => Try[WomValue]): ValueEvaluator[A] =
    forOperationWithShortCircuit(op, PartialFunction.empty[WomValue, WomValue])


  private def forOperationWithShortCircuit[A <: BinaryOperation](op: (WomValue, WomValue) => Try[WomValue],
                                                                 shortCircuit: PartialFunction[WomValue, WomValue]) = new ValueEvaluator[A] {
    override def evaluateValue(a: A,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      a.left.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions) flatMap { left =>
        if (shortCircuit.isDefinedAt(left.value)) {
          EvaluatedValue(shortCircuit(left.value), left.sideEffectFiles).validNel
        } else {
          a.right.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions) flatMap { right =>
            val rawResult = op(left.value, right.value)

            // Allow unsupplied optionals, but only if we're instantiating a command:
            val handleOptionals = rawResult.recover {
              case OptionalNotSuppliedException(_) if forCommandInstantiationOptions.isDefined => WomString("")
            }

            handleOptionals.toErrorOr map { newValue => EvaluatedValue(newValue, left.sideEffectFiles ++ right.sideEffectFiles) }
          }
        }
      }
  }
}
