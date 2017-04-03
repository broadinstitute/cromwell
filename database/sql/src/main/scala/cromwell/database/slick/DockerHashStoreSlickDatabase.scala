package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.DockerHashStoreSqlDatabase
import cromwell.database.sql.tables.DockerHashStoreEntry

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

trait DockerHashStoreSlickDatabase extends DockerHashStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  /**
    * Adds a docker hash entry to the store.
    */
  override def addDockerHashStoreEntry(dockerHashStoreEntry: DockerHashStoreEntry)
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.dockerHashStoreEntries += dockerHashStoreEntry
    runTransaction(action) void
  }

  /**
    * Retrieves docker hash entries for a workflow.
    *
    */
  override def queryDockerHashStoreEntries(workflowExecutionUuid: String)
                                          (implicit ec: ExecutionContext): Future[Seq[DockerHashStoreEntry]] = {
    val action = dataAccess.dockerHashStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  /**
    * Deletes docker hash entries related to a workflow, returning the number of rows affected.
    */
  override def removeDockerHashStoreEntries(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.dockerHashStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }
}
