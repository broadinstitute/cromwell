package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.types.WomType

@typeclass
trait TypeEvaluator[A] {
  def evaluateType(a: A,
                   linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                  (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType]
}
