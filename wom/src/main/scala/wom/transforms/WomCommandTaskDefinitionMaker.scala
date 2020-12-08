package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.CommandTaskDefinition
import simulacrum._

@typeclass
trait WomCommandTaskDefinitionMaker[A] {
  def toWomTaskDefinition(a: A): ErrorOr[CommandTaskDefinition]
}
