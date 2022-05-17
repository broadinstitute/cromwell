package centaur.test.workflow

import cats.effect.IO
import cats.instances.list._
import cats.syntax.traverse._
import cromwell.api.model.{SubmittedWorkflow, WorkflowId}

/**
 * Tracks submitted workflow ids to enable cleanup should a test need to be retried.
 */
class SubmittedWorkflowTracker {
  private var submittedWorkflowIds: List[WorkflowId] = List.empty

  /**
   * Run the specified cleanup function on the submitted workflow IDs tracked by this `CentaurTestCase`, clearing out
   * the list of submitted workflow IDs afterward.
   */
  def cleanUpBeforeRetry(cleanUpFunction: WorkflowId => IO[Unit]): IO[Unit] = for {
    _ <- submittedWorkflowIds.traverse(cleanUpFunction)
    _ = submittedWorkflowIds = List.empty
  } yield ()

  /**
   * Add a `SubmittedWorkflow` to the list of `SubmittedWorkflow`s to clean up should the test case represented by this
   * object require a retry. Prevents unwanted cache hits from partially successful attempts when retrying a call
   * caching test case.
   */
  def add(submittedWorkflow: SubmittedWorkflow): Unit = {
    submittedWorkflowIds = submittedWorkflow.id :: submittedWorkflowIds
  }
}
