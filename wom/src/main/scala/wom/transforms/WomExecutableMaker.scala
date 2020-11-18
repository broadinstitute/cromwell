package wom.transforms

import common.Checked
import wom.executable.{Executable, WomBundle}
import simulacrum._
import wom.core.WorkflowJson
import wom.expression.IoFunctionSet

@typeclass
trait WomExecutableMaker[A] {
  def toWomExecutable(a: A, inputs: Option[WorkflowJson], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable]
}

@typeclass
trait WomBundleMaker[A] {
  def toWomBundle(a: A): Checked[WomBundle]
}
