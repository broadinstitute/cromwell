package wdl.model.draft3.graph

import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wom.callable.Callable
import wom.types.WomType

import scala.language.implicitConversions

@typeclass
trait GraphElementValueConsumer[A] {
  def graphElementConsumedValueHooks(a: A, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[UnlinkedConsumedValueHook]]
}

@typeclass
trait ExpressionValueConsumer[A] {
  def expressionConsumedValueHooks(a: A): Set[UnlinkedConsumedValueHook]
}
