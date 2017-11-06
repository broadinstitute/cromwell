package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def updateWorkflowsInState(updates: List[(String, String)])
                                     (implicit ec: ExecutionContext): Future[Unit] = {
    val updateQueries = updates.map({
      case (oldState, newState) => dataAccess.workflowStateForWorkflowState(oldState).update(newState)
    })

    val action = DBIO.sequence(updateQueries)
    runTransaction(action) void
  }

  override def updateToRestartedIfInState(states: List[String])
                                         (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      _ <- dataAccess.workflowRestartedIfStateIsIn(states).update(true)
    } yield ()

    runTransaction(action) void
  }

  override def updateWorkflowState(workflowId: String, newWorkflowState: String)
                                     (implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    val action =  for {
      restarted <- dataAccess.workflowRestartedForId(workflowId).result.headOption
      _ <- dataAccess.workflowStateForId(workflowId).update(newWorkflowState)
    } yield restarted
    
    runTransaction(action)
  }

  override def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreEntryIdsAutoInc ++= workflowStoreEntries
    runTransaction(action) void
  }

  override def queryWorkflowStoreEntries(limit: Int, queryWorkflowState: String, queryWorkflowRestarted: Boolean, updateWorkflowState: String)
                                        (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]] = {
    val action = for {
      workflowStoreEntries <- dataAccess.workflowStoreEntriesForWorkflowState((queryWorkflowState, queryWorkflowRestarted, limit.toLong)).result
      _ <- DBIO.sequence(workflowStoreEntries map updateWorkflowStateAndRestartedForWorkflowExecutionUuid(updateWorkflowState, restarted = false))
    } yield workflowStoreEntries
    runTransaction(action)
  }

  private def updateWorkflowStateAndRestartedForWorkflowExecutionUuid(updateWorkflowState: String, restarted: Boolean)
                                                         (workflowStoreEntry: WorkflowStoreEntry)
                                                         (implicit ec: ExecutionContext): DBIO[Unit] = {
    val workflowExecutionUuid = workflowStoreEntry.workflowExecutionUuid
    for {
      updateCount <- dataAccess.workflowStateAndRestartedForWorkflowExecutionUuid(workflowExecutionUuid).update((updateWorkflowState, restarted))
      _ <- assertUpdateCount(s"Update $workflowExecutionUuid to $updateWorkflowState", updateCount, 1)
    } yield ()
  }

  override def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }
  
  override def stats(implicit ec: ExecutionContext): Future[Map[String, Int]] = {
    val action = dataAccess.workflowStoreStats.result
    runTransaction(action) map { _.toMap }
  }
}
