package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.WomExpression
import wom.types.WomType

@typeclass
trait WomExpressionMaker[A] {
  def makeWomExpression(a: A,
                        typeAliases: Map[String, WomType],
                        consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomExpression]
}
