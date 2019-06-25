package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import cromwell.core.{HogGroup, WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreAbortResponse.WorkflowStoreAbortResponse
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowStoreAbortResponse, WorkflowStoreState, WorkflowSubmissionResponse}

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

class InMemoryWorkflowStore extends WorkflowStore {

  var workflowStore = Map.empty[WorkflowIdAndSources, WorkflowStoreState]

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowSubmissionResponse]] = {
    val actualWorkflowState = if (sources.head.workflowOnHold) WorkflowStoreState.OnHold else WorkflowStoreState.Submitted
    val addedWorkflows = sources map { WorkflowIdAndSources(WorkflowId.randomId(), _) -> actualWorkflowState }
    workflowStore ++= addedWorkflows.toList.toMap
    Future.successful(addedWorkflows map {
      case (WorkflowIdAndSources(id, _), _) => WorkflowSubmissionResponse(actualWorkflowState, id)
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
      case (workflow, WorkflowStoreState.Submitted) =>
        val hogGroup = HogGroup.decide(workflow.sources.workflowOptions, workflow.id)
        WorkflowToStart(workflow.id, OffsetDateTime.now, workflow.sources, Submitted, hogGroup)
      case _ => throw new IllegalArgumentException("This workflow is not currently in a startable state")
    }

    Future.successful(workflowsToStart.toList)
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

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStoreAbortResponse] = {
    workflowStore collectFirst {
      case (workflowIdAndSources, workflowStoreState) if workflowIdAndSources.id == id =>
        (workflowIdAndSources, workflowStoreState)
    } match {
      case Some((workflowIdAndSources, WorkflowStoreState.OnHold)) =>
        workflowStore -= workflowIdAndSources
        Future.successful(WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted)
      case Some((workflowIdAndSources, _)) =>
        workflowStore += workflowIdAndSources -> WorkflowStoreState.Aborting
        // In memory workflows can never be restarted (since this is destroyed on a server restart)
        Future.successful(WorkflowStoreAbortResponse.AbortRequested)
      case None =>
        Future.successful(WorkflowStoreAbortResponse.NotFound)
    }
  }

  override def writeWorkflowHeartbeats(workflowIds: Set[(WorkflowId, OffsetDateTime)],
                                       heartbeatDateTime: OffsetDateTime)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    Future.successful(workflowIds.size)
  }

  override def switchOnHoldToSubmitted(id: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())

  override def findWorkflowsWithAbortRequested(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[WorkflowId]] = Future.successful(List.empty)

  override def findWorkflows(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[WorkflowId]] = Future.successful(workflowStore.keys.map(_.id))
}

final case class WorkflowIdAndSources(id: WorkflowId, sources: WorkflowSourceFilesCollection)
