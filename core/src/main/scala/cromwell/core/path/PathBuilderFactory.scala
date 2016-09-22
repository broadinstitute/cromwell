package cromwell.core.path

import cromwell.core.WorkflowOptions

/**
  * Provide a method that can instantiate a path builder with the specified workflow options.
  */
trait PathBuilderFactory {
  /**
    * Typically this should be called once per workflow.
    */
  def withOptions(options: WorkflowOptions): PathBuilder
}
