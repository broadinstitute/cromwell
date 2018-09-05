package cromwell.database.sql

import cromwell.database.sql.tables.JobKeyValueEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobKeyValueSqlDatabase {
  this: SqlDatabase =>

  def addJobKeyValueEntry(jobKeyValueEntry: JobKeyValueEntry)
                         (implicit ec: ExecutionContext): Future[Unit]

  def addJobKeyValueEntries(jobKeyValueEntries: Iterable[JobKeyValueEntry])
                           (implicit ec: ExecutionContext): Future[Unit]

  def queryJobKeyValueEntries(workflowExecutionUuid: String)
                             (implicit ec: ExecutionContext): Future[Seq[JobKeyValueEntry]]

  def queryStoreValue(workflowExecutionUuid: String, callFqn: String, jobScatterIndex: Int,
                      jobRetryAttempt: Int, storeKey: String)
                     (implicit ec: ExecutionContext): Future[Option[String]]
}
