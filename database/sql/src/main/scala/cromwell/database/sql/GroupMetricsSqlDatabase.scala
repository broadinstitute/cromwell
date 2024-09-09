package cromwell.database.sql

import cromwell.database.sql.tables.GroupMetricsEntry

import java.sql.Timestamp
import scala.concurrent.{ExecutionContext, Future}

trait GroupMetricsSqlDatabase {

  this: SqlDatabase =>

  /**
   * Insert or update Group Metrics entry to the table
   *
   */
  def recordGroupMetricsEntry(groupMetricsEntry: GroupMetricsEntry)(implicit ec: ExecutionContext): Future[Unit]

  /**
   * Returns number of entries associated with given group
   */
  def countGroupMetricsEntries(groupId: String)(implicit ec: ExecutionContext): Future[Int]

  /**
   * Returns list of groups **currently** experiencing quota exhaustion
   */
  def getQuotaExhaustedGroups(thresholdTimestamp: Timestamp)(implicit ec: ExecutionContext): Future[Seq[String]]
}
