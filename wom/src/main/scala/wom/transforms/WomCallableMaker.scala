package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomCallableMaker[A] {
  def toWomCallable(a: A): ErrorOr[Callable]
}
