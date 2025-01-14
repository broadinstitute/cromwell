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

    /*
      The approach here is to try and insert the record into the table and if that fails because a record with that
      group_id already exists, then it will update that record with new quota_exhaustion_detected timestamp.
      The insert and update happen in 2 different transactions. Using 2 separate transactions avoids deadlocks that
      occur during upserts or when using Slick's insertOrUpdate in a single transaction. Deadlocks occur in scenarios
      with concurrent threads trying to update the table, even under the stricter Serializable transaction isolation
      level. This is caused by a gap lock on the IX_GROUP_METRICS_ENTRY_GI index, which locks gaps between index
      records on a page, leading to deadlocks for transactions involving the same or nearby group_id values.
      See https://broadworkbench.atlassian.net/browse/AN-286 and https://broadworkbench.atlassian.net/browse/WX-1847.

      This approach should have minimal performance impact since the table has only 3 columns and low cardinality.

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
