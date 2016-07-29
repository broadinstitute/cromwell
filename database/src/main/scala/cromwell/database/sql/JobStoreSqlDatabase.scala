package cromwell.database.sql

import cromwell.database.sql.tables.JobStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreSqlDatabase {
  this: SqlDatabase =>

  def addJobStoreEntries(entries: Iterable[JobStoreEntry])(implicit ec: ExecutionContext): Future[Unit]

  def fetchJobResult(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                            (implicit ec: ExecutionContext): Future[Option[JobStoreEntry]]

  def removeJobResultsForWorkflows(workflowUuids: Iterable[String])(implicit ec: ExecutionContext): Future[Int]
}
