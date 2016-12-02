package cromwell.core

import wdl4s.{WdlJson, WdlSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */

sealed trait WorkflowSourceFilesCollection {
  def wdlSource: WdlSource
  def inputsJson: WdlJson
  def workflowOptionsJson: WorkflowOptionsJson
  def importsZipFileOption: Option[Array[Byte]] = this match {
    case _: WorkflowSourceFilesWithoutImports => None
    case WorkflowSourceFilesWithDependenciesZip(_, _, _, importsZip) => Option(importsZip) // i.e. Some(importsZip) if our wiring is correct
  }

  def copyOptions(workflowOptions: WorkflowOptionsJson) = this match {
    case w: WorkflowSourceFilesWithoutImports => WorkflowSourceFilesWithoutImports(w.wdlSource, w.inputsJson, workflowOptions)
    case w: WorkflowSourceFilesWithDependenciesZip => WorkflowSourceFilesWithDependenciesZip(w.wdlSource, w.inputsJson, workflowOptions, w.importsZip)
  }
}

object WorkflowSourceFilesCollection {
  def apply(wdlSource: WdlSource,
            inputsJson: WdlJson,
            workflowOptionsJson: WorkflowOptionsJson,
            importsFile: Option[Array[Byte]]): WorkflowSourceFilesCollection = importsFile match {
    case Some(imports) => WorkflowSourceFilesWithDependenciesZip(wdlSource, inputsJson, workflowOptionsJson, imports)
    case None => WorkflowSourceFilesWithoutImports(wdlSource, inputsJson, workflowOptionsJson)
  }
}

final case class WorkflowSourceFilesWithoutImports(wdlSource: WdlSource,
                                                   inputsJson: WdlJson,
                                                   workflowOptionsJson: WorkflowOptionsJson) extends WorkflowSourceFilesCollection

final case class WorkflowSourceFilesWithDependenciesZip(wdlSource: WdlSource,
                                                        inputsJson: WdlJson,
                                                        workflowOptionsJson: WorkflowOptionsJson,
                                                        importsZip: Array[Byte]) extends WorkflowSourceFilesCollection {
  override def toString = s"WorkflowSourceFilesWithDependenciesZip($wdlSource, $inputsJson, $workflowOptionsJson, <<ZIP BINARY CONTENT>>)"
}
