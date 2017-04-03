package cromwell.database.sql

import cromwell.database.sql.tables.DockerHashStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait DockerHashStoreSqlDatabase {
  this: SqlDatabase =>

  /**
    * Adds a docker hash entry to the store.
    *
    */
  def addDockerHashStoreEntry(dockerHashStoreEntry: DockerHashStoreEntry)
                             (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Retrieves docker hash entries for a workflow.
    *
    */
  def queryDockerHashStoreEntries(workflowExecutionUuid: String)
                                 (implicit ec: ExecutionContext): Future[Seq[DockerHashStoreEntry]]

  /**
    * Deletes docker hash entries related to a workflow, returning the number of rows affected.
    */
  def removeDockerHashStoreEntries(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]
}
