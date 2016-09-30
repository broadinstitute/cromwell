package cromwell.database.slick

import cromwell.database.sql.tables.SummaryStatusEntry

import scala.concurrent.ExecutionContext

trait SummaryStatusSlickDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def getSummaryStatusEntryMaximumId(summaryTableName: String, summarizedTableName: String)
                                                   (implicit ec: ExecutionContext): DBIO[Option[Long]] = {
    dataAccess.
      maximumIdForSummaryTableNameSummarizedTableName((summaryTableName, summarizedTableName)).
      result.headOption
  }

  private[slick] def previousOrMaximum(previous: Long, longs: Seq[Long]): Long = (previous +: longs).max

  private[slick] def upsertSummaryStatusEntryMaximumId(summaryTableName: String, summarizedTableName: String,
                                                       maximumId: Long)(implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.summaryStatusEntryIdsAutoInc.
          insertOrUpdate(SummaryStatusEntry(summaryTableName, summarizedTableName, maximumId))
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          maximumIdForSummaryTableNameSummarizedTableName((summaryTableName, summarizedTableName)).
          update(maximumId)
        _ <- updateCount match {
          case 0 =>
            dataAccess.summaryStatusEntryIdsAutoInc +=
              SummaryStatusEntry(summaryTableName, summarizedTableName, maximumId)
          case _ => assertUpdateCount("upsertSummaryStatusEntryMaximumId", updateCount, 1)
        }
      } yield ()
    }
  }
}
