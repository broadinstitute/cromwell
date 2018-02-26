package wdl.model.draft3.graph

import simulacrum.typeclass

import scala.language.implicitConversions

@typeclass
trait UnlinkedValueConsumer[A] {
  def consumedValueHooks(a: A): Set[UnlinkedConsumedValueHook]
}
