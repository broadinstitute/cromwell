package cromwell.database.sql

import cromwell.database.sql.tables.SubWorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait SubWorkflowStoreSqlDatabase {
  this: SqlDatabase =>

  def addSubWorkflowStoreEntry(rootWorkflowExecutionUuid: String,
                               parentWorkflowExecutionUuid: String,
                               callFullyQualifiedName: String,
                               jobIndex: Int,
                               jobAttempt: Int,
                               subWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Unit]

  def querySubWorkflowStore(parentWorkflowExecutionUuid: String, callFqn: String, jobIndex: Int, jobAttempt: Int)
                    (implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]]

  def removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]
}
