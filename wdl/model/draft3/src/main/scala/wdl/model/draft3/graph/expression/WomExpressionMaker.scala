package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr

import scala.language.implicitConversions
import simulacrum.typeclass
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.WomExpression

@typeclass
trait WomExpressionMaker[A] {
  def makeWomExpression(a: A,
                        consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomExpression]
}
