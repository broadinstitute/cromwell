package cromwell.backend

import cromwell.backend.model.{WorkflowDescriptor, JobDescriptor}

/**
  * Defines basic structure and functionality to work with a backend.
  */
trait Backend {
  /**
    * Defines needed data to be able to execute a workflow.
    */
  val workflowDescriptor: WorkflowDescriptor

  /**
    * Registers code to be executed before the backend is ready for executing jobs for the specific workflow.
    */
  def beforeAll(): Unit

  /**
    * Factory for creating Jobs in the specific backend.
    *
    * @param jobDescriptor All information needed to execute a task in the backend.
    * @return A BackendJobExecutor with functionality to handle the life cycle management of the task.
    */
  def createJobExecutor(jobDescriptor: JobDescriptor): BackendJobExecutor

  /**
    * Registers code to be executed after the backend finished executing all related tasks for the specific workflow.
    */
  def afterAll(): Unit

  /**
    * Executes validation on workflow descriptor in order to see if the workflow can be executed by the backend.
    */
  def validate(): Unit
}
