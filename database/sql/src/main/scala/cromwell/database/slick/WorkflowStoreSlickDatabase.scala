package cromwell.database.slick

import java.sql.Timestamp
import java.text.DateFormat

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def setAllRunningToAborting()(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess
      .workflowStateForWorkflowState(WorkflowStoreState.Running)
      .update(WorkflowStoreState.Aborting)

    runTransaction(action).void
  }

  override def markRunningAndAbortingAsRestarted(cromwellId: Option[String])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.heartbeatTimestamp(cromwellId).update(None)

    runTransaction(action).void
  }

  /**
    * Set the workflow Id to Aborting state.
    * @return Some(restarted) if the workflow exists in the store, where restarted is true if the workflow's heartbeat
    *         timestamp was cleared on restart.
    *         None if the workflow does not exist in the store
    */
  override def setToAborting(workflowId: String)
                            (implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    val action =  for {
      restarted <- dataAccess.heartbeatClearedForWorkflowId(workflowId).result.headOption
      _ <- dataAccess.workflowStateForId(workflowId).update(WorkflowStoreState.Aborting)
    } yield restarted
    
    runTransaction(action)
  }

  override def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreEntryIdsAutoInc ++= workflowStoreEntries
    runTransaction(action).void
  }

  private def now: Timestamp = new Timestamp(System.currentTimeMillis())

  override def fetchStartableWorkflows(limit: Int, cromwellId: Option[String], heartbeatTtl: FiniteDuration)
                                      (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]] = {
    val action = for {
      workflowStoreEntries <- dataAccess.fetchStartableWorkflows((limit.toLong, cromwellId, heartbeatTtl.ago)).result
      _ <- DBIO.sequence(workflowStoreEntries map updateForFetched(cromwellId))
    } yield workflowStoreEntries

    runTransaction(action)
  }

  private def updateForFetched(cromwellId: Option[String])(workflowStoreEntry: WorkflowStoreEntry)
                                                          (implicit ec: ExecutionContext): DBIO[Unit] = {
    val workflowExecutionUuid = workflowStoreEntry.workflowExecutionUuid
    val updateState = workflowStoreEntry.workflowState match {
        // Submitted workflows become running when fetched
      case WorkflowStoreState.Submitted => WorkflowStoreState.Running
        // Running or Aborting stay as is
      case other => other
    }

    val displayNow = DateFormat.getDateTimeInstance.format(now)

    for {
      // When fetched, the heartbeat timestamp is set to now so we don't pick it up next time.
      updateCount <- dataAccess.workflowStoreFieldsForPickup(workflowExecutionUuid).update((updateState, cromwellId, Option(now)))
      _ <- assertUpdateCount(s"Update $workflowExecutionUuid to $updateState, heartbeat timestamp to $displayNow", updateCount, 1)
    } yield ()
  }

  override def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }
  
  override def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]] = {
    val action = dataAccess.workflowStoreStats.result
    runTransaction(action) map { _.toMap }
  }
}
