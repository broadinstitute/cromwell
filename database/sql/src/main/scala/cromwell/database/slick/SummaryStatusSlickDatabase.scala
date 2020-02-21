package cromwell.database.slick

import cromwell.database.sql.tables.SummaryStatusEntry

import scala.concurrent.ExecutionContext

trait SummaryStatusSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def getSummaryStatusEntrySummaryPosition(summaryName: String): DBIO[Option[Long]] = {
    dataAccess.summaryPositionForSummaryName(summaryName).result.headOption
  }

  private[slick] def upsertSummaryStatusEntrySummaryPosition(summaryName: String,
                                                             summaryPosition: Long)
                                                            (implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.summaryStatusEntryIdsAutoInc.
          insertOrUpdate(SummaryStatusEntry(summaryName, summaryPosition))
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.summaryPositionForSummaryName(summaryName).update(summaryPosition)
        _ <- updateCount match {
          case 0 =>
            dataAccess.summaryStatusEntryIdsAutoInc +=
              SummaryStatusEntry(summaryName, summaryPosition)
          case _ => assertUpdateCount("upsertSummaryStatusEntrySummaryPosition", updateCount, 1)
        }
      } yield ()
    }
  }
}
