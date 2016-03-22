package cromwell.backend

import cromwell.backend.model.JobDescriptor

/**
  * Defines basic structure and functionality to execute a job in a backend.
  */
trait BackendJobExecutor {
  /**
    * Defines needed data to be able to execute a job.
    */
  val jobDescriptor: JobDescriptor

  /**
    * Prepare the task and context for execution.
    */
  def prepare(): Unit

  /**
    * Executes task in given context.
    */
  def execute(): Unit

  /**
    * Stops a task execution.
    */
  def stop(): Unit

  /**
    * Performs a cleanUp after the task was executed.
    */
  def cleanUp(): Unit
}
