package wdl.draft3.transforms.wdlom2wom.linking

import simulacrum.typeclass
import scala.language.implicitConversions

@typeclass
trait UnlinkedValueConsumer[A] {
  def consumedValueNames(a: A): Set[UnlinkedConsumedValueName]
}

