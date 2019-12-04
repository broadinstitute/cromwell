package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.sql.tables.SummaryQueueEntry

import scala.concurrent.ExecutionContext

trait SummaryQueueSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def writeSummaryQueueEntry(metadataEntryIds: Seq[Long], timestamp: Timestamp) = {
    dataAccess.summaryQueueEntryIdsAutoInc ++= metadataEntryIds.map(id => SummaryQueueEntry(id, timestamp, None, None))
  }

}
