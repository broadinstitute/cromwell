package cromwell.core

import wdl4s.{WdlJson, WdlSource}

/**
  * Represents the collection of source files that a user submits to run a workflow
  */
final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson,
                                     workflowOptionsJson: WorkflowOptionsJson)
