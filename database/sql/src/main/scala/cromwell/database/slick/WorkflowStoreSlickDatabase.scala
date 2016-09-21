package cromwell.database.slick

import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def updateWorkflowState(queryWorkflowState: String, updateWorkflowState: String)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStateForWorkflowState(queryWorkflowState).update(updateWorkflowState)
    runTransaction(action) map { _ => () }
  }

  override def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreEntryIdsAutoInc ++= workflowStoreEntries
    runTransaction(action) map { _ => () }
  }

  override def queryWorkflowStoreEntries(limit: Int, queryWorkflowState: String, updateWorkflowState: String)
                                        (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]] = {
    val action = for {
      workflowStoreEntries <- dataAccess.workflowStoreEntriesForWorkflowState(queryWorkflowState, limit).result
      _ <- DBIO.sequence(workflowStoreEntries map updateWorkflowStateForWorkflowExecutionUuid(updateWorkflowState))
    } yield workflowStoreEntries
    runTransaction(action)
  }

  private def updateWorkflowStateForWorkflowExecutionUuid(updateWorkflowState: String)
                                                         (workflowStoreEntry: WorkflowStoreEntry)
                                                         (implicit ec: ExecutionContext): DBIO[Unit] = {
    val workflowExecutionUuid = workflowStoreEntry.workflowExecutionUuid
    for {
      updateCount <- dataAccess.workflowStateForWorkflowExecutionUuid(workflowExecutionUuid).update(updateWorkflowState)
      _ <- assertUpdateCount(s"Update $workflowExecutionUuid to $updateWorkflowState", updateCount, 1)
    } yield ()
  }

  override def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }
}
