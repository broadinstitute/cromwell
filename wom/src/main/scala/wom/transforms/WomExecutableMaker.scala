package wom.transforms

import common.Checked
import wom.executable.Executable

trait WomExecutableMaker[A] {
  def toWomExecutable(a: A, inputFile: Option[String]): Checked[Executable]
}

object WomExecutableMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomExecutableMaker[A]): WomExecutableMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeExecutable[A: WomExecutableMaker](val a: A) {
    def toWomExecutable(inputs: Option[String]) = WomExecutableMaker[A].toWomExecutable(a, inputs)
    def toWomExecutable = WomExecutableMaker[A].toWomExecutable(a, None)
  }
}
