package wdl.transforms.draft2.wdlom2wom

import common.Checked
import wdl.WdlNamespaceWithWorkflow
import wom.executable.Executable
import wom.transforms.WomExecutableMaker

object WdlDraft2WomExecutableMaker extends WomExecutableMaker[WdlNamespaceWithWorkflow] {
  override def toWomExecutable(wdlNamespaceWithWorkflow: WdlNamespaceWithWorkflow, inputFile: Option[String] = None): Checked[Executable] = {
    WdlDraft2InputParsing.buildWomExecutable(wdlNamespaceWithWorkflow.workflow, inputFile)
  }
}
