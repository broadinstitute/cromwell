package cromwell.subworkflowstore

import cromwell.database.sql.tables.SubWorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait SubWorkflowStore {
  def addSubWorkflowStoreEntry(subWorkflowStoreEntry: SubWorkflowStoreEntry)(implicit ec: ExecutionContext): Future[Unit]

  def querySubWorkflowStore(parentWorkflowExecutionUuid: String, callFqn: String, jobIndex: Int, jobAttempt: Int)
                           (implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]]

  def removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]
}
