package cromwell.engine.workflow.workflowstore

import cats.data.NonEmptyList
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowSubmissionResponse

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

class InMemoryWorkflowStore extends WorkflowStore {

  var workflowStore = Map.empty[SubmittedWorkflow, WorkflowStoreState]

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowSubmissionResponse]] = {
    val actualWorkflowState = if (sources.head.workflowOnHold) WorkflowStoreState.OnHold else WorkflowStoreState.Submitted
    val submittedWorkflows = sources map { SubmittedWorkflow(WorkflowId.randomId(), _) -> actualWorkflowState }
    workflowStore = workflowStore ++ submittedWorkflows.toList.toMap
    Future.successful(submittedWorkflows map {
      case (SubmittedWorkflow(id, _), _) => WorkflowSubmissionResponse(actualWorkflowState, id)
    })
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    val startableWorkflows = workflowStore filter { _._2 == WorkflowStoreState.Submitted } take n
    val updatedWorkflows = startableWorkflows map { _._1 -> WorkflowStoreState.Running }
    workflowStore = workflowStore ++ updatedWorkflows

    val workflowsToStart = startableWorkflows map {
      case (workflow, WorkflowStoreState.Submitted) => WorkflowToStart(workflow.id, workflow.sources, Submitted)
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

  override def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]] = Future.successful(Map(WorkflowStoreState.Submitted -> workflowStore.size))

  override def abortAllRunning()(implicit ec: ExecutionContext): Future[Unit] = {
    workflowStore = workflowStore.map({
      case (workflow, WorkflowStoreState.Running) => workflow -> WorkflowStoreState.Aborting
      case (workflow, state) => workflow -> state
    })
    Future.successful(())
  }

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    if (workflowStore.exists(_._1.id == id)) {
      val state = workflowStore.find(_._1.id == id)
      workflowStore = workflowStore ++ workflowStore.find(_._1.id == id).map({ _._1 -> WorkflowStoreState.Aborting }).toMap
      // In memory workflows can never be restarted (since this is destroyed on a server restart)
      Future.successful(state.map(_ => false))
    } else {
      Future.successful(None)
    }
  }

  override def writeWorkflowHeartbeats(workflowIds: Set[WorkflowId])(implicit ec: ExecutionContext): Future[Int] =
    Future.successful(workflowIds.size)

  override def switchOnHoldToSubmitted(id: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())
}

final case class SubmittedWorkflow(id: WorkflowId, sources: WorkflowSourceFilesCollection)
