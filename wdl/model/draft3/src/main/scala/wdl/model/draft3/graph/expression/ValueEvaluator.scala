package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.IoFunctionSet
import wom.values.WomValue

import scala.language.implicitConversions

@typeclass
trait ValueEvaluator[A] {
  def evaluateValue(a: A,
                    inputs: Map[String, WomValue],
                    ioFunctionSet: IoFunctionSet,
                    linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomValue]
}
