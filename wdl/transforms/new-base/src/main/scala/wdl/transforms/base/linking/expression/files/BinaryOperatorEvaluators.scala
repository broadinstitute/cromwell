package wdl.transforms.base.linking.expression.files

import cats.syntax.apply._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}

import scala.util.Try

object BinaryOperatorEvaluators {
  implicit val logicalOrEvaluator: FileEvaluator[LogicalOr] = forOperation(_.or(_))
  implicit val logicalAndEvaluator: FileEvaluator[LogicalAnd] = forOperation(_.and(_))
  implicit val equalsEvaluator: FileEvaluator[Equals] = forOperation(_.equals(_))
  implicit val notEqualsEvaluator: FileEvaluator[NotEquals] = forOperation(_.notEquals(_))
  implicit val lessThanEvaluator: FileEvaluator[LessThan] = forOperation(_.lessThan(_))
  implicit val lessThanOrEqualEvaluator: FileEvaluator[LessThanOrEquals] = forOperation(_.lessThanOrEqual(_))
  implicit val greaterThanEvaluator: FileEvaluator[GreaterThan] = forOperation(_.greaterThan(_))
  implicit val greaterThanOrEqualEvaluator: FileEvaluator[GreaterThanOrEquals] = forOperation(_.greaterThanOrEqual(_))
  implicit val addEvaluator: FileEvaluator[Add] = forOperation(_.add(_))
  implicit val subtractEvaluator: FileEvaluator[Subtract] = forOperation(_.subtract(_))
  implicit val multiplyEvaluator: FileEvaluator[Multiply] = forOperation(_.multiply(_))
  implicit val divideEvaluator: FileEvaluator[Divide] = forOperation(_.divide(_))
  implicit val remainderEvaluator: FileEvaluator[Remainder] = forOperation(_.mod(_))

  private def forOperation[A <: BinaryOperation](op: (WomValue, WomValue) => Try[WomValue]) = new FileEvaluator[A] {
    override def predictFilesNeededToEvaluate(a: A,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType
    )(implicit
      fileEvaluator: FileEvaluator[ExpressionElement],
      valueEvaluator: ValueEvaluator[ExpressionElement]
    ): ErrorOr[Set[WomFile]] =
      (a.left.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
       a.right.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
      ) mapN { _ ++ _ }
  }
}
