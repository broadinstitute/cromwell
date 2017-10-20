package cromwell.subworkflowstore
import cromwell.database.sql.SubWorkflowStoreSqlDatabase
import cromwell.database.sql.tables.SubWorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

class SqlSubWorkflowStore(subWorkflowStoreSqlDatabase: SubWorkflowStoreSqlDatabase) extends SubWorkflowStore {
  override def addSubWorkflowStoreEntry(rootWorkflowExecutionUuid: String,
                                        parentWorkflowExecutionUuid: String,
                                        callFullyQualifiedName: String,
                                        jobIndex: Int,
                                        jobAttempt: Int,
                                        subWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Unit] = {
    subWorkflowStoreSqlDatabase.addSubWorkflowStoreEntry(
      rootWorkflowExecutionUuid,
      parentWorkflowExecutionUuid,
      callFullyQualifiedName,
      jobIndex,
      jobAttempt,
      subWorkflowExecutionUuid
    )
  }

  override def querySubWorkflowStore(parentWorkflowExecutionUuid: String, callFqn: String, jobIndex: Int, jobAttempt: Int)(implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]] = {
    subWorkflowStoreSqlDatabase.querySubWorkflowStore(parentWorkflowExecutionUuid, callFqn, jobIndex, jobAttempt)
  }

  override def removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    subWorkflowStoreSqlDatabase.removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid)
  }
}
