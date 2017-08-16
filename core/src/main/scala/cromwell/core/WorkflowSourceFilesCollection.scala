package cromwell.core

import wdl4s.wdl.{WorkflowJson, WorkflowSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */

sealed trait WorkflowSourceFilesCollection {
  def workflowSource: WorkflowSource
  def inputsJson: WorkflowJson
  def workflowOptionsJson: WorkflowOptionsJson
  def labelsJson: WorkflowJson
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
  def apply(workflowSource: WorkflowSource,
            workflowType: Option[WorkflowType],
            workflowTypeVersion: Option[WorkflowTypeVersion],
            inputsJson: WorkflowJson,
            workflowOptionsJson: WorkflowOptionsJson,
            labelsJson: WorkflowJson,
            importsFile: Option[Array[Byte]]): WorkflowSourceFilesCollection = importsFile match {
    case Some(imports) =>
      WorkflowSourceFilesWithDependenciesZip(workflowSource, workflowType, workflowTypeVersion, inputsJson, workflowOptionsJson, labelsJson, imports)
    case None =>
      WorkflowSourceFilesWithoutImports(workflowSource, workflowType, workflowTypeVersion, inputsJson, workflowOptionsJson, labelsJson)
  }
}

final case class WorkflowSourceFilesWithoutImports(workflowSource: WorkflowSource,
                                                   workflowType: Option[WorkflowType],
                                                   workflowTypeVersion: Option[WorkflowTypeVersion],
                                                   inputsJson: WorkflowJson,
                                                   workflowOptionsJson: WorkflowOptionsJson,
                                                   labelsJson: WorkflowJson) extends WorkflowSourceFilesCollection

final case class WorkflowSourceFilesWithDependenciesZip(workflowSource: WorkflowSource,
                                                        workflowType: Option[WorkflowType],
                                                        workflowTypeVersion: Option[WorkflowTypeVersion],
                                                        inputsJson: WorkflowJson,
                                                        workflowOptionsJson: WorkflowOptionsJson,
                                                        labelsJson: WorkflowJson,
                                                        importsZip: Array[Byte]) extends WorkflowSourceFilesCollection {
  override def toString = s"WorkflowSourceFilesWithDependenciesZip($workflowSource, $inputsJson, $workflowOptionsJson, $labelsJson, <<ZIP BINARY CONTENT>>)"
}
