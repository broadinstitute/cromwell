package wdl.model.draft3.graph.expression

import common.validation.ErrorOr.ErrorOr

import scala.language.implicitConversions
import simulacrum.typeclass
import wom.types.WomType

@typeclass
trait WomTypeMaker[A] {
  def determineWomType(a: A, availableAliases: Map[String, WomType]): ErrorOr[WomType]
}
