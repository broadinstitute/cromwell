package cromwell.database.slick

import cromwell.database.sql.tables.SummaryStatusEntry

import scala.concurrent.ExecutionContext

trait SummaryStatusSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def getSummaryStatusEntrySummaryPosition(summaryName: String): DBIO[Option[Long]] = {
    dataAccess.summaryPositionForSummaryName(summaryName).result.headOption
  }

  /* If you're about to (re-)introduce slick upserts (or insertOrUpdates):
    * See https://github.com/broadinstitute/cromwell/pull/5332 which removed them for failing tests in Slick 3.3.2
    * Make sure the new slick version you're using passes KeyValueDatabaseSpec (and all the others, of course)
 */
  private[slick] def upsertSummaryStatusEntrySummaryPosition(summaryName: String,
                                                             summaryPosition: Long)
                                                            (implicit ec: ExecutionContext): DBIO[Unit] = {
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
