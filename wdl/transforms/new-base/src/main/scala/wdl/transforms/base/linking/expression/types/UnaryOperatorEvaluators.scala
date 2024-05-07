package wdl.transforms.base.linking.expression.types

import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types.WomType

import scala.util.Try

object UnaryOperatorEvaluators {

  implicit val unaryPlusEvaluator: TypeEvaluator[UnaryPlus] = forOperation(_.unaryPlus)
  implicit val unaryNegationEvaluator: TypeEvaluator[UnaryNegation] = forOperation(_.unaryMinus)
  implicit val logicalNotEvaluator: TypeEvaluator[LogicalNot] = forOperation(_.not)

  private def forOperation[A <: UnaryOperation](op: WomType => Try[WomType]) = new TypeEvaluator[A] {
    override def evaluateType(a: A,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      a.argument.evaluateType(linkedValues, typeAliases) flatMap { op(_).toErrorOr }
  }
}
