package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.SubWorkflowStoreSqlDatabase
import cromwell.database.sql.tables.SubWorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

trait SubWorkflowStoreSlickDatabase extends SubWorkflowStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  def addSubWorkflowStoreEntry(subWorkflowStoreEntry: SubWorkflowStoreEntry)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      subWorkflowStoreEntryId <- dataAccess.subWorkflowStoreEntryIdsAutoInc += subWorkflowStoreEntry
    } yield ()
    runTransaction(action) void
  }

  override def querySubWorkflowStore(parentWorkflowExecutionUuid: String, callFqn: String, jobIndex: Int, jobAttempt: Int)
                                    (implicit ec: ExecutionContext): Future[Option[SubWorkflowStoreEntry]] = {
    val action = for {
      subWorkflowStoreEntryOption <- dataAccess.subWorkflowStoreEntriesForJobKey(
        (parentWorkflowExecutionUuid, callFqn, jobIndex, jobAttempt)
      ).result.headOption
    } yield subWorkflowStoreEntryOption

    runTransaction(action)
  }

  override def removeSubWorkflowStoreEntries(parentWorkflowExecutionUuid: String)
                                            (implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      deleted <- dataAccess.subWorkflowStoreEntriesForParentWorkflowExecutionUuid(parentWorkflowExecutionUuid).delete
    } yield deleted
    
    runTransaction(action)
  }
}
