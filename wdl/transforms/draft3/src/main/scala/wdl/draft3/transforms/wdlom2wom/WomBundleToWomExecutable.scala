package wdl.draft3.transforms.wdlom2wom

import common.Checked
import wdl.shared.transforms.wdlom2wom.WdlSharedInputParsing
import wom.core.WorkflowJson
import wom.executable.{Executable, WomBundle}
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker

object WomBundleToWomExecutable {

  implicit val draft3WomBundleToWomExecutable: WomExecutableMaker[WomBundle] = new WomExecutableMaker[WomBundle] {
    override def toWomExecutable(a: WomBundle, inputs: Option[WorkflowJson], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = WdlSharedInputParsing.buildWomExecutable(a, inputs, ioFunctions, strictValidation)
  }
}
