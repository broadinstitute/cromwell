package cromwell.core

import better.files.File
import wdl4s.{WdlJson, WdlSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */

sealed trait WorkflowSourceFilesCollection {
  def wdlSource: WdlSource
  def inputsJson: WdlJson
  def workflowOptionsJson: WorkflowOptionsJson

  def copyOptions(workflowOptions: WorkflowOptionsJson) = this match {
    case w: WorkflowSourceFiles => WorkflowSourceFiles(w.wdlSource, w.inputsJson, workflowOptions)
    case w: WorkflowSourceFilesWithImports => WorkflowSourceFilesWithImports(w.wdlSource, w.inputsJson, workflowOptions, w.importsFile)
  }

}

final case class WorkflowSourceFiles(wdlSource: WdlSource,
                                     inputsJson: WdlJson,
                                     workflowOptionsJson: WorkflowOptionsJson) extends WorkflowSourceFilesCollection

final case class WorkflowSourceFilesWithImports(wdlSource: WdlSource,
                                                inputsJson: WdlJson,
                                                workflowOptionsJson: WorkflowOptionsJson,
                                                importsFile: File) extends WorkflowSourceFilesCollection
