package cromwell.engine.workflow.workflowstore

import cats.data.NonEmptyList
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

import scala.concurrent.{ExecutionContext, Future}

class InMemoryWorkflowStore extends WorkflowStore {

  var workflowStore = List.empty[SubmittedWorkflow]

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {
    val submittedWorkflows = sources map { SubmittedWorkflow(WorkflowId.randomId(), _, WorkflowStoreState.Submitted) }
    workflowStore = workflowStore ++ submittedWorkflows.toList
    Future.successful(submittedWorkflows map { _.id })
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    val startableWorkflows = workflowStore filter { _.state == state } take n
    val updatedWorkflows = startableWorkflows map { _.copy(state = WorkflowStoreState.Running) }
    workflowStore = (workflowStore diff startableWorkflows) ++ updatedWorkflows

    Future.successful(startableWorkflows map { _.toWorkflowToStart })
  }

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    if (workflowStore.exists(_.id == id)) {
      workflowStore = workflowStore filterNot { _.id == id }
      Future.successful(true)
    } else {
      Future.successful(false)
    }
  }

  override def initialize(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())
}

final case class SubmittedWorkflow(id: WorkflowId, sources: WorkflowSourceFilesCollection, state: WorkflowStoreState) {
  def toWorkflowToStart: WorkflowToStart = {
    state match {
      case r: StartableState => WorkflowToStart(id, sources, r)
      case _ => throw new IllegalArgumentException("This workflow is not currently in a startable state")
    }
  }
}