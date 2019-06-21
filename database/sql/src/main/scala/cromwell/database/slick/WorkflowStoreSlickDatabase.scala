package cromwell.database.slick

import java.sql.Timestamp

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import slick.jdbc.TransactionIsolation

import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSlickDatabase extends WorkflowStoreSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def setStateToState(fromWorkflowState: String, toWorkflowState: String)
                              (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess
      .workflowStateForWorkflowState(fromWorkflowState)
      .update(toWorkflowState)

    runTransaction(action).void
  }

  override def deleteOrUpdateWorkflowToState(workflowExecutionUuid: String,
                                             workflowStateToDelete1: String,
                                             workflowStateToDelete2: String,
                                             workflowStateForUpdate: String)
                                            (implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    val action = for {
      // First, delete all rows in either of our states to be deleted.
      deleted <- dataAccess
        .workflowStoreEntryForWorkflowExecutionUUidAndWorkflowStates(
          (workflowExecutionUuid, workflowStateToDelete1, workflowStateToDelete2)
        )
        .delete
      // Second, for any single row still present update its state.
      updated <-
        deleted match {
          case 0 =>
            dataAccess.workflowStateForWorkflowExecutionUUid(workflowExecutionUuid).update(workflowStateForUpdate)
          case _ => assertUpdateCount("deleteOrUpdateWorkflowToState", deleted, 1)
        }
    } yield (deleted, updated)

    runTransaction(action) map { case (deleted, updated) => if (deleted == 0 && updated == 0) None else Option(deleted > 0) }
  }

  override def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.workflowStoreEntryIdsAutoInc ++= workflowStoreEntries
    runTransaction(action).void
  }

  override def fetchWorkflowsInState(limit: Int,
                                     cromwellId: String,
                                     heartbeatTimestampTimedOut: Timestamp,
                                     heartbeatTimestampTo: Timestamp,
                                     workflowStateFrom: String,
                                     workflowStateTo: String,
                                     workflowStateExcluded: String)
                                    (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]] = {

    def updateForFetched(cromwellId: String,
                         heartbeatTimestampTo: Timestamp,
                         workflowStateFrom: String,
                         workflowStateTo: String)
                        (workflowStoreEntry: WorkflowStoreEntry)
                        (implicit ec: ExecutionContext): DBIO[Unit] = {
      val workflowExecutionUuid = workflowStoreEntry.workflowExecutionUuid
      val updateState = workflowStoreEntry.workflowState match {
        case matched if matched == workflowStateFrom => workflowStateTo
        case other => other
      }

      for {
        // When fetched, the heartbeat timestamp is set to now so we don't pick it up next time.
        updateCount <- dataAccess
          .workflowStoreFieldsForPickup(workflowExecutionUuid)
          .update((updateState, Option(cromwellId), Option(heartbeatTimestampTo)))
        _ <- assertUpdateCount(
          s"Update $workflowExecutionUuid to $updateState, heartbeat timestamp to $heartbeatTimestampTo",
          updateCount,
          1
        )
      } yield ()
    }

    val action = for {
      workflowStoreEntries <- dataAccess.fetchStartableWorkflows(
        (limit.toLong, heartbeatTimestampTimedOut, workflowStateExcluded)
      ).result
      _ <- DBIO.sequence(
        workflowStoreEntries map updateForFetched(cromwellId, heartbeatTimestampTo, workflowStateFrom, workflowStateTo)
      )
    } yield workflowStoreEntries

    // This should be safe because there are no repeated reads and the same rows are being locked for update as with
    // the default RepeatableRead isolation. This has the advantage of avoiding MySQL's aggressive Serializable-ish
    // gap and next-key locking behavior in RepeatableRead mode with non-index scans so now only the rows that are
    // modified will be locked. Any rows that are inserted due to the absence of the gap and next-key locks should be
    // valid candidates for workflow pickup.
    // https://dev.mysql.com/doc/refman/8.0/en/innodb-transaction-isolation-levels.html
    runTransaction(action, TransactionIsolation.ReadCommitted)
  }

  override def writeWorkflowHeartbeats(workflowExecutionUuids: Seq[String],
                                       heartbeatTimestamp: Timestamp)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    // Return the count of heartbeats written. This could legitimately be less than the size of the `workflowExecutionUuids`
    // List if any of those workflows completed and their workflow store entries were removed.
    val action = for {
      counts <- DBIO.sequence(workflowExecutionUuids map { workflowExecutionUuid =>
        dataAccess.heartbeatForWorkflowStoreEntry(workflowExecutionUuid).update(Option(heartbeatTimestamp))
      })
    } yield counts.sum
    // Auto-commit mode, so each statement is individually committed to avoid deadlocks
    runAction(action)
  }

  override def releaseWorkflowStoreEntries(cromwellId: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.releaseWorkflowStoreEntries(cromwellId).update((None, None))
    runTransaction(action)
  }

  override def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.workflowStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    runTransaction(action)
  }

  override def workflowStateCounts(implicit ec: ExecutionContext): Future[Map[String, Int]] = {
    val action = dataAccess.workflowStoreStats.result
    runTransaction(action) map { _.toMap }
  }

  override def updateWorkflowState(workflowExecutionUuid: String,
                                   fromWorkflowState: String,
                                   toWorkflowState: String)
                                  (implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      updated <- dataAccess
        .workflowStateForWorkflowExecutionUUidAndWorkflowState((workflowExecutionUuid, fromWorkflowState))
        .update(toWorkflowState)
    } yield updated

    runTransaction(action)
  }

  override def findWorkflowsWithAbortRequested(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[String]] = {
    runTransaction(dataAccess.findWorkflowsWithAbortRequested(cromwellId).result)
  }

  override def findWorkflows(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[String]] = {
    runTransaction(dataAccess.findWorkflows(cromwellId).result)
  }
}
