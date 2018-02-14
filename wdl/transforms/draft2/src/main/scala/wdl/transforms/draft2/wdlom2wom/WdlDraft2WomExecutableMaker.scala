package wdl.transforms.draft2.wdlom2wom

import common.Checked
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wom.executable.Executable
import wom.transforms.WomExecutableMaker
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

object WdlDraft2WomExecutableMaker extends WomExecutableMaker[WdlNamespaceWithWorkflow] {
  override def toWomExecutable(inputs: ExecutableMakerInputs[WdlNamespaceWithWorkflow]): Checked[Executable] = {
    WdlDraft2InputParsing.buildWomExecutable(inputs.from.workflow, inputs.inputs)
  }

  override def toWomExecutable(a: WdlNamespaceWithWorkflow, inputs: Option[String]): Checked[Executable] = toWomExecutable(ExecutableMakerInputs(a, inputs))
}
