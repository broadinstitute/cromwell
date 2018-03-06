package wdl.draft3.transforms.linking.expression.values

import common.validation.Validation._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.IoFunctionSet
import wom.values.WomValue

import scala.util.Try

object UnaryOperatorEvaluators {

  implicit val unaryPlusEvaluator: ValueEvaluator[UnaryPlus] = forOperation(_.unaryPlus)
  implicit val unaryNegationEvaluator: ValueEvaluator[UnaryNegation] = forOperation(_.unaryMinus)
  implicit val logicalNotEvaluator: ValueEvaluator[LogicalNot] = forOperation(_.not)

  private def forOperation[A <: UnaryOperation](op: WomValue => Try[WomValue]) = new ValueEvaluator[A] {
    override def evaluateValue(a: A,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      a.argument.evaluateValue(inputs, ioFunctionSet) flatMap { op(_).toErrorOr }
  }
}
