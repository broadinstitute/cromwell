package centaur.test

import centaur.test.workflow.Workflow
import cromwell.api.model.{SubmittedWorkflow, WorkflowMetadata}

/**
  * An exception with information about a centaur test failure.
  *
  * @param message A messaged returned for the exception.
  * @param testName The ScalaTest name that was being executed.
  * @param workflowIdOption The optional workflow id.
  * @param metadataJsonOption The optional metadata.
  * @param causeOption The optional underlying cause.
  */
case class CentaurTestException private(message: String,
                                        testName: String,
                                        workflowIdOption: Option[String],
                                        metadataJsonOption: Option[String],
                                        causeOption: Option[Exception])
  extends RuntimeException(message, causeOption.orNull)

object CentaurTestException {

  /** Create a new CentaurTestException for a completed workflow. */
  def apply(message: String,
            workflowDefinition: Workflow,
            submittedWorkflow: SubmittedWorkflow,
            actualMetadata: WorkflowMetadata): CentaurTestException = {
    new CentaurTestException(
      message,
      workflowDefinition.testName,
      Option(submittedWorkflow.id.toString),
      Option(actualMetadata.value),
      None
    )
  }

  /** Create a new CentaurTestException for a submitted workflow. */
  def apply(message: String,
            workflowDefinition: Workflow,
            submittedWorkflow: SubmittedWorkflow): CentaurTestException = {
    new CentaurTestException(
      message,
      workflowDefinition.testName,
      Option(submittedWorkflow.id.toString),
      None,
      None
    )
  }

  /** Create a new CentaurTestException for only a workflow definition. */
  def apply(message: String, workflowDefinition: Workflow): CentaurTestException = {
    new CentaurTestException(
      message,
      workflowDefinition.testName,
      None,
      None,
      None
    )
  }

  /** Create a new CentaurTestException for only a workflow definition, including a root cause. */
  def apply(message: String, workflowDefinition: Workflow, cause: Exception): CentaurTestException = {
    new CentaurTestException(
      message,
      workflowDefinition.testName,
      None,
      None,
      Option(cause)
    )
  }
}
