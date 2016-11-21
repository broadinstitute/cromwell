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

  def addSubWorkflowStoreEntry(rootWorkflowExecutionUuid: String, 
                               parentWorkflowExecutionUuid: String,
                               callFullyQualifiedName: String,
                               jobIndex: Int,
                               jobAttempt: Int,
                               subWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowStoreEntry <- dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(rootWorkflowExecutionUuid).result.headOption
      _ <- workflowStoreEntry match {
        case Some(rootWorkflow) =>
          dataAccess.subWorkflowStoreEntryIdsAutoInc +=
            SubWorkflowStoreEntry(
              rootWorkflow.workflowStoreEntryId,
              parentWorkflowExecutionUuid,
              callFullyQualifiedName,
              jobIndex,
              jobAttempt,
              subWorkflowExecutionUuid
            )
        case None => DBIO.failed(new IllegalArgumentException(s"Could not find root workflow with UUID $rootWorkflowExecutionUuid"))
      }
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

  override def removeSubWorkflowStoreEntries(rootWorkflowExecutionUuid: String)
                                            (implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      workflowStoreEntry <- dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(rootWorkflowExecutionUuid).result.headOption
      deleted <- workflowStoreEntry match {
        case Some(rootWorkflow) =>
          dataAccess.subWorkflowStoreEntriesForRootWorkflowId(rootWorkflow.workflowStoreEntryId.get).delete
        case None =>
          DBIO.successful(0)
      }
    } yield deleted
    
    runTransaction(action)
  }
}
