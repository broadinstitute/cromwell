package cromwell.database.slick

import java.sql.Timestamp
import java.time.OffsetDateTime
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.SummaryQueueEntry

trait SummaryQueueSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def writeSummaryQueueEntry(metadataEntryIds: Seq[Long], timestamp: Timestamp) = {
    dataAccess.summaryQueueEntryIdsAutoInc ++= metadataEntryIds.map(id => SummaryQueueEntry(id, timestamp, None, None))
  }

  private[slick] def markSummaryQueueEntriesAsSummarizedForMetadataEntryIds(metadataEntryIds: Seq[Long]) = {
    val updateTimestamp = OffsetDateTime.now.toSystemTimestamp

    DBIO.sequence(metadataEntryIds map { id =>
      dataAccess.summaryQueueSummarizationTimestampForId(id).update(Option(updateTimestamp))
    })
  }



}
