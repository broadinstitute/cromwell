package cromwell.database.slick

import cromwell.database.sql.GroupMetricsSqlDatabase
import cromwell.database.sql.tables.GroupMetricsEntry

import scala.concurrent.{ExecutionContext, Future}

trait GroupMetricsSlickDatabase extends GroupMetricsSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def recordGroupMetricsEntry(
    groupMetricsEntry: GroupMetricsEntry
  )(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      updateCount <- dataAccess
        .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
        .update(groupMetricsEntry.quotaExhaustionDetected)
      _ <- updateCount match {
        case 0 => dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
        case _ => assertUpdateCount("recordGroupMetricsEntry", updateCount, 1)
      }
    } yield ()
    runTransaction(action)
  }

  override def countGroupMetricsEntries(groupId: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countGroupMetricsEntriesForGroupId(groupId).result
    runTransaction(action)
  }
}
