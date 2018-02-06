package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable

trait WomCallableMaker[A] {
  def toWomCallable(a: A): ErrorOr[Callable]
}

object WomCallableMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomCallableMaker[A]): WomCallableMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeCallable[A: WomCallableMaker](val a: A) {
    def toWomCallable = WomCallableMaker[A].toWomCallable(a)
  }
}
