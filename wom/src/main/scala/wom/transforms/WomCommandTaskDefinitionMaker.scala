package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.CommandTaskDefinition

trait WomCommandTaskDefinitionMaker[A] {
  def toWomTaskDefinition(a: A): ErrorOr[CommandTaskDefinition]
}

object WomCommandTaskDefinitionMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomCommandTaskDefinitionMaker[A]): WomCommandTaskDefinitionMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeTaskDefinition[A: WomCommandTaskDefinitionMaker](val a: A) {
    def toWomTaskDefinition = WomCommandTaskDefinitionMaker[A].toWomTaskDefinition(a)
  }
}
