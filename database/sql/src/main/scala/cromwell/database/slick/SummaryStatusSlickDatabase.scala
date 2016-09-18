package cromwell.database.slick

import cromwell.database.sql.tables.SummaryStatus

import scala.concurrent.ExecutionContext

trait SummarizingSlickDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def getSummaryStatusMaximumId(summaryTableName: String, summarizedTableName: String)
                                              (implicit ec: ExecutionContext): DBIO[Option[Long]] = {
    dataAccess.
      summaryStatusMaximumIdBySummaryTableNameSummarizedTableName(summaryTableName, summarizedTableName).
      result.headOption
  }

  private[slick] def maximumOrZero(longs: Seq[Long]): Long = (0L +: longs).max

  private[slick] def upsertSummaryStatusMaximumId(summaryTableName: String, summarizedTableName: String,
                                                  maximumId: Long)(implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.summaryStatusIdsAutoInc.
          insertOrUpdate(SummaryStatus(summaryTableName, summarizedTableName, maximumId))
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          summaryStatusMaximumIdBySummaryTableNameSummarizedTableName(summaryTableName, summarizedTableName).
          update(maximumId)
        _ <- updateCount match {
          case 0 =>
            dataAccess.summaryStatusIdsAutoInc += SummaryStatus(summaryTableName, summarizedTableName, maximumId)
          case _ => assertUpdateCount("upsertSummaryStatusMaximumId", updateCount, 1)
        }
      } yield ()
    }
  }
}
