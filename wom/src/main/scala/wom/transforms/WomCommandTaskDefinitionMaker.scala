package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.CommandTaskDefinition
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomCommandTaskDefinitionMaker[A] {
  @op("toWomTaskDefinition")
  def toWomTaskDefinition(a: A): ErrorOr[CommandTaskDefinition]
}
