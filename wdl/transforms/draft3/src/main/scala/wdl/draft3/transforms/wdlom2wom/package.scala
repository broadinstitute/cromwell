package wdl.draft3.transforms

import better.files.File
import wdl.draft3.transforms.ast2wdlom._
import wom.executable.Executable

package object wdlom2wom {
  implicit val fromWorkflowDefinitionElementToWorkflowDefinition = WorkflowDefinitionElementToWomWorkflowDefinition

  def womFromDraft3FileConverter(inputs: Option[String]): FromAtoB[File, Executable] = FromAtoB.viaX(draft3FileElementFromFile, FileElementToWomExecutable(inputs))
}
