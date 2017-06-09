package cromwell.core

import wdl4s.{WdlJson, WdlSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */

sealed trait WorkflowSourceFilesCollection {
  def workflowSource: WdlSource
  def inputsJson: WdlJson
  def workflowOptionsJson: WorkflowOptionsJson
  def labelsJson: WdlJson
  def workflowType: Option[WorkflowType]
  def workflowTypeVersion: Option[WorkflowTypeVersion]

  def importsZipFileOption: Option[Array[Byte]] = this match {
    case _: WorkflowSourceFilesWithoutImports => None
    case w: WorkflowSourceFilesWithDependenciesZip => Option(w.importsZip) // i.e. Some(importsZip) if our wiring is correct
  }

  def copyOptions(workflowOptions: WorkflowOptionsJson) = this match {
    case w: WorkflowSourceFilesWithoutImports => w.copy(workflowOptionsJson = workflowOptions)
    case w: WorkflowSourceFilesWithDependenciesZip => w.copy(workflowOptionsJson = workflowOptions)
  }
}

object WorkflowSourceFilesCollection {
  def apply(workflowSource: WdlSource,
            workflowType: Option[WorkflowType],
            workflowTypeVersion: Option[WorkflowTypeVersion],
            inputsJson: WdlJson,
            workflowOptionsJson: WorkflowOptionsJson,
            labelsJson: WdlJson,
            importsFile: Option[Array[Byte]]): WorkflowSourceFilesCollection = importsFile match {
    case Some(imports) =>
      WorkflowSourceFilesWithDependenciesZip(workflowSource, workflowType, workflowTypeVersion, inputsJson, workflowOptionsJson, labelsJson, imports)
    case None =>
      WorkflowSourceFilesWithoutImports(workflowSource, workflowType, workflowTypeVersion, inputsJson, workflowOptionsJson, labelsJson)
  }
}

final case class WorkflowSourceFilesWithoutImports(workflowSource: WdlSource,
                                                   workflowType: Option[WorkflowType],
                                                   workflowTypeVersion: Option[WorkflowTypeVersion],
                                                   inputsJson: WdlJson,
                                                   workflowOptionsJson: WorkflowOptionsJson,
                                                   labelsJson: WdlJson) extends WorkflowSourceFilesCollection

final case class WorkflowSourceFilesWithDependenciesZip(workflowSource: WdlSource,
                                                        workflowType: Option[WorkflowType],
                                                        workflowTypeVersion: Option[WorkflowTypeVersion],
                                                        inputsJson: WdlJson,
                                                        workflowOptionsJson: WorkflowOptionsJson,
                                                        labelsJson: WdlJson,
                                                        importsZip: Array[Byte]) extends WorkflowSourceFilesCollection {
  override def toString = s"WorkflowSourceFilesWithDependenciesZip($workflowSource, $inputsJson, $workflowOptionsJson, $labelsJson, <<ZIP BINARY CONTENT>>)"
}
