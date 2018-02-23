package wdl.draft3.transforms.wdlom2wom.linking

import scala.language.implicitConversions

import simulacrum.typeclass

@typeclass
trait UnlinkedValueGenerator[A] {
  def generatedValueNames(a: A): Set[UnlinkedGeneratedValueName]
}

