package cromwell.subworkflowstore
import cromwell.database.sql.SubWorkflowStoreSqlDatabase
import cromwell.database.sql.tables.SubWorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

class SqlSubWorkflowStore(subWorkflowStoreSqlDatabase: SubWorkflowStoreSqlDatabase) extends SubWorkflowStore {
  override def addSubWorkflowStoreEntry(subWorkflowStoreEntry: SubWorkflowStoreEntry)(implicit ec: ExecutionContext): Future[Unit] = {
    subWorkflowStoreSqlDatabase.addSubWorkflowStoreEntry(subWorkflowStoreEntry)
  }

  override def querySubWorkflowStore(parentWorkflowExecutionUuid: String, callFqn: String, jobIndex: Int, jobAttempt: Int)(implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]] = {
    subWorkflowStoreSqlDatabase.querySubWorkflowStore(parentWorkflowExecutionUuid, callFqn, jobIndex, jobAttempt)
  }

  override def removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    subWorkflowStoreSqlDatabase.removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid)
  }
}
