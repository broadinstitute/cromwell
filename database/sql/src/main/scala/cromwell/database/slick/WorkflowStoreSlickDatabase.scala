package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def setAllRunningToAborting()(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess
      .workflowStateForWorkflowState(WorkflowStoreState.Running)
      .update(WorkflowStoreState.Aborting)

    runTransaction(action).void
  }

  override def markRunningAndAbortingAsRestarted()(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.restartedFlagForRunningAndAborting.update(true)

    runTransaction(action).void
  }

  /**
    * Set the workflow Id to Aborting state.
    * @return Some(restarted) if the workflow exists in the store, where restarted is the value of its restarted flag
    *         None if the workflow does not exist in the store
    */
  override def setToAborting(workflowId: String)
                            (implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    val action =  for {
      restarted <- dataAccess.workflowRestartedForId(workflowId).result.headOption
      _ <- dataAccess.workflowStateForId(workflowId).update(WorkflowStoreState.Aborting)
    } yield restarted
    
    runTransaction(action)
  }

  override def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreEntryIdsAutoInc ++= workflowStoreEntries
    runTransaction(action).void
  }

  override def fetchStartableWorkflows(limit: Int)
                                      (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]] = {
    val action = for {
      workflowStoreEntries <- dataAccess.fetchStartableWorkflows(limit.toLong).result
      _ <- DBIO.sequence(workflowStoreEntries map updateWorkflowStateAndRestartedForWorkflowExecutionUuid)
    } yield workflowStoreEntries

    runTransaction(action)
  }

  private def updateWorkflowStateAndRestartedForWorkflowExecutionUuid(workflowStoreEntry: WorkflowStoreEntry)
                                                         (implicit ec: ExecutionContext): DBIO[Unit] = {
    val workflowExecutionUuid = workflowStoreEntry.workflowExecutionUuid
    val updateState = workflowStoreEntry.workflowState match {
        // Submitted workflows become running when fetched
      case WorkflowStoreState.Submitted => WorkflowStoreState.Running
        // Running or Aborting stay as is
      case other => other
    }
    for {
      // When fetched, the restarted flag is set back to false so we don't pick it up next time.
      updateCount <- dataAccess.workflowStateAndRestartedForWorkflowExecutionUuid(workflowExecutionUuid).update((updateState, false))
      _ <- assertUpdateCount(s"Update $workflowExecutionUuid to $updateState", updateCount, 1)
    } yield ()
  }

  override def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }
  
  override def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]] = {
    val action = dataAccess.workflowStoreStats.result
    runTransaction(action) map { _.toMap }
  }
}
