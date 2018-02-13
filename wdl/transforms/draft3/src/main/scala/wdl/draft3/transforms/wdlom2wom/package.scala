package wdl.draft3.transforms

import better.files.File
import cats.instances.either._

import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.wdlom2wom.FileElementToWomExecutable.FileElementAndInputsFile
import wom.executable.Executable

package object wdlom2wom {
  implicit val fromWorkflowDefinitionElementToWorkflowDefinition = WorkflowDefinitionElementToWomWorkflowDefinition

  def womFromDraft3FileConverter(inputs: Option[String]): CheckedAtoB[File, Executable] = {
    draft3FileElementFromFile.map(FileElementAndInputsFile(_, inputs)) andThen FileElementToWomExecutable.instance
  }
}
