package cromwell.database.slick

import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def updateWorkflowState(oldState: String, newState: String)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreStateByWorkflowStoreState(oldState).update(newState)
    runTransaction(action) map { _ => () }
  }

  override def addWorkflow(entries: Iterable[WorkflowStoreEntry])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreAutoInc ++= entries
    runTransaction(action) map { _ => () }
  }

  override def fetchRunnableWorkflows(limit: Int, fetchState: String, updateState: String)
                                     (implicit ec: ExecutionContext):
  Future[Seq[WorkflowStoreEntry]] = {

    val action = for {
      workflowEntries <- dataAccess.workflowStoreEntriesByState(fetchState, limit).result
      _ <- DBIO.sequence(workflowEntries map verifyUpdate(updateState))
    } yield workflowEntries
    runTransaction(action)
  }

  private def verifyUpdate(updateState: String)
                          (workflowStoreEntry: WorkflowStoreEntry)
                          (implicit ec: ExecutionContext): DBIO[Unit] = {
    val workflowUuid = workflowStoreEntry.workflowUuid
    for {
      updateCount <- dataAccess.workflowStateByUuid(workflowUuid).update(updateState)
      _ <- assertUpdateCount(s"Update $workflowUuid to $updateState", updateCount, 1)
    } yield ()
  }

  override def removeWorkflow(id: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntryByWorkflowUuid(id).delete
    runTransaction(action)
  }
}
