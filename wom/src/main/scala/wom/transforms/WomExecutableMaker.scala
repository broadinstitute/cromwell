package wom.transforms

import common.Checked
import wom.executable.Executable
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomExecutableMaker[A] {
  def toWomExecutable(a: A, inputFile: Option[String] = None): Checked[Executable]
}
