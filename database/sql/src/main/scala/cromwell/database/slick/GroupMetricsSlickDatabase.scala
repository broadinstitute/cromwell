package cromwell.database.slick

import cromwell.database.sql.GroupMetricsSqlDatabase
import cromwell.database.sql.tables.GroupMetricsEntry
import org.postgresql.util.PSQLException
import slick.jdbc.TransactionIsolation

import java.sql.{SQLIntegrityConstraintViolationException, Timestamp}
import scala.concurrent.{ExecutionContext, Future}

trait GroupMetricsSlickDatabase extends GroupMetricsSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def recordGroupMetricsEntry(
    groupMetricsEntry: GroupMetricsEntry
  )(implicit ec: ExecutionContext): Future[Unit] = {
    def updateGroupMetricsEntry(): Future[Unit] = {
      val updateAction = for {
        _ <- dataAccess
          .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
          .update(groupMetricsEntry.quotaExhaustionDetected)
      } yield ()
      runTransaction(updateAction, TransactionIsolation.ReadCommitted)
    }

    /* The approach here is to try and insert the record into database in 1 transaction and if that fails, because
       a record with that group_id already exists, then it will update that record with new quota_exhaustion_detected
       timestamp in a separate transaction. The reason for 2 separate transactions are because when manual upsert or
       Slick's insertOrUpdate is performed in a single transaction, and there are threads trying to update the table
       concurrently it results in a deadlock (even with stricter transaction isolation level Serializable). This
       happens because it gets a gap lock on index IX_GROUP_METRICS_ENTRY_GI which locks the gaps between index records
       within the page and as a result it runs into a deadlock when there are transactions trying to insert/update for
       same group ID or for group IDs that exist in the locked records.
       See https://broadworkbench.atlassian.net/browse/AN-286 and https://broadworkbench.atlassian.net/browse/WX-1847.

       Note: a unique constraint on group_id also exists in database.
     */
    val insertAction = for {
      _ <- dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
    } yield ()

    runTransaction(insertAction, TransactionIsolation.ReadCommitted)
      .recoverWith {
        case ex
            if ex.isInstanceOf[SQLIntegrityConstraintViolationException] || (ex
              .isInstanceOf[PSQLException] && ex.getMessage.contains(
              "duplicate key value violates unique constraint \"UC_GROUP_METRICS_ENTRY_GI\""
            )) =>
          updateGroupMetricsEntry()
        case e => Future.failed(e)
      }
  }

  override def countGroupMetricsEntries(groupId: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countGroupMetricsEntriesForGroupId(groupId).result
    runTransaction(action)
  }

  override def getQuotaExhaustedGroups(
    thresholdTimestamp: Timestamp
  )(implicit ec: ExecutionContext): Future[Seq[String]] = {
    val action = dataAccess.groupsExperiencingQuotaExhaustion(thresholdTimestamp).result
    runTransaction(action)
  }
}
