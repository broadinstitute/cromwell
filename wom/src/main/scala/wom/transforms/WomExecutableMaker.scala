package wom.transforms

import common.Checked
import wom.executable.Executable
import simulacrum._
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

import scala.language.implicitConversions

@typeclass
trait WomExecutableMaker[A] {
  def toWomExecutable(inputs: ExecutableMakerInputs[A]): Checked[Executable]
  def toWomExecutable(a: A, inputs: Option[String] = None): Checked[Executable]
}

object WomExecutableMaker {
  final case class ExecutableMakerInputs[A](from: A, inputs: Option[String])
}