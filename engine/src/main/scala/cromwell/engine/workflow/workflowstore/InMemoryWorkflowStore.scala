package cromwell.engine.workflow.workflowstore

import cats.data.NonEmptyList
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.{Aborting, Running, StartableState}

import scala.concurrent.{ExecutionContext, Future}

class InMemoryWorkflowStore extends WorkflowStore {

  var workflowStore = Map.empty[SubmittedWorkflow, WorkflowStoreState]

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {
    val submittedWorkflows = sources map { SubmittedWorkflow(WorkflowId.randomId(), _) -> WorkflowStoreState.Submitted }
    workflowStore = workflowStore ++ submittedWorkflows.toList.toMap
    Future.successful(submittedWorkflows map { _._1.id })
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    val startableWorkflows = workflowStore filter { _._2 == state } take n
    val updatedWorkflows = startableWorkflows map { _._1 -> WorkflowStoreState.Running }
    workflowStore = workflowStore ++ updatedWorkflows

    val workflowsToStart = startableWorkflows map {
      case (workflow, s: StartableState) => WorkflowToStart(workflow.id, workflow.sources, s)
      case _ => throw new IllegalArgumentException("This workflow is not currently in a startable state")
    }

    Future.successful(workflowsToStart.toList)
  }

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    if (workflowStore.exists(_._1.id == id)) {
      workflowStore = workflowStore filterNot { _._1.id == id }
      Future.successful(true)
    } else {
      Future.successful(false)
    }
  }

  override def initialize(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())

  override def stats(implicit ec: ExecutionContext): Future[Map[String, Int]] = Future.successful(Map("Submitted" -> workflowStore.size))

  override def abortAllRunning()(implicit ec: ExecutionContext): Future[Unit] = {
    workflowStore = workflowStore.map({
      case (workflow, Running) => workflow -> Aborting
      case (workflow, state) => workflow -> state
    })
    Future.successful(())
  }
  
  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    if (workflowStore.exists(_._1.id == id)) {
      workflowStore = workflowStore ++ workflowStore.find(_._1.id == id).map({ _._1 -> Aborting }).toMap
      Future.successful(true)
    } else {
      Future.successful(false)
    }
  }
}

final case class SubmittedWorkflow(id: WorkflowId, sources: WorkflowSourceFilesCollection)
