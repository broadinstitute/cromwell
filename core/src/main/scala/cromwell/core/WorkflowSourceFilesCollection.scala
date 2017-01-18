package cromwell.core

import wdl4s.{WdlJson, WdlSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */

sealed trait WorkflowSourceFilesCollection {
  def wdlSource: WdlSource
  def inputsJson: WdlJson
  def workflowOptionsJson: WorkflowOptionsJson
  def labelsJson: WdlJson


  def importsZipFileOption: Option[Array[Byte]] = this match {
    case _: WorkflowSourceFilesWithoutImports => None
    case WorkflowSourceFilesWithDependenciesZip(_, _, _, _, importsZip) => Option(importsZip) // i.e. Some(importsZip) if our wiring is correct
  }

  def copyOptions(workflowOptions: WorkflowOptionsJson) = this match {
    case w: WorkflowSourceFilesWithoutImports => WorkflowSourceFilesWithoutImports(
      wdlSource = w.wdlSource,
      inputsJson = w.inputsJson,
      workflowOptionsJson = workflowOptions,
      labelsJson = w.labelsJson)

    case w: WorkflowSourceFilesWithDependenciesZip => WorkflowSourceFilesWithDependenciesZip(
      wdlSource = w.wdlSource,
      inputsJson = w.inputsJson,
      workflowOptionsJson = workflowOptions,
      labelsJson = w.labelsJson,
      importsZip = w.importsZip)
  }
}

object WorkflowSourceFilesCollection {
  def apply(wdlSource: WdlSource,
            inputsJson: WdlJson,
            workflowOptionsJson: WorkflowOptionsJson,
            labelsJson: WdlJson,
            importsFile: Option[Array[Byte]]): WorkflowSourceFilesCollection = importsFile match {
    case Some(imports) => WorkflowSourceFilesWithDependenciesZip(wdlSource, inputsJson, workflowOptionsJson, labelsJson, imports)
    case None => WorkflowSourceFilesWithoutImports(wdlSource, inputsJson, workflowOptionsJson, labelsJson)
  }
}

final case class WorkflowSourceFilesWithoutImports(wdlSource: WdlSource,
                                                   inputsJson: WdlJson,
                                                   workflowOptionsJson: WorkflowOptionsJson,
                                                   labelsJson: WdlJson) extends WorkflowSourceFilesCollection

final case class WorkflowSourceFilesWithDependenciesZip(wdlSource: WdlSource,
                                                        inputsJson: WdlJson,
                                                        workflowOptionsJson: WorkflowOptionsJson,
                                                        labelsJson: WdlJson,
                                                        importsZip: Array[Byte]) extends WorkflowSourceFilesCollection {
  override def toString = s"WorkflowSourceFilesWithDependenciesZip($wdlSource, $inputsJson, $workflowOptionsJson, $labelsJson, <<ZIP BINARY CONTENT>>)"
}
