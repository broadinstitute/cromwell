package wdl.model.draft3.graph

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wom.types.WomType

import scala.language.implicitConversions

@typeclass
trait GraphElementValueConsumer[A] {
  def graphElementConsumedValueHooks(a: A, typeAliases: Map[String, WomType]): ErrorOr[Set[UnlinkedConsumedValueHook]]
}

@typeclass
trait ExpressionValueConsumer[A] {
  def expressionConsumedValueHooks(a: A): Set[UnlinkedConsumedValueHook]
}
