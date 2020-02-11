package cromwell.database.slick

import cromwell.database.sql.tables.SummaryQueueEntry

trait SummaryQueueSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def writeSummaryQueueEntries(metadataJournalIds: Seq[Long]) = {
    dataAccess.summaryQueueEntries ++= metadataJournalIds.map(id => SummaryQueueEntry(id))
  }

  private[slick] def deleteSummaryQueueEntriesByMetadataJournalIds(metadataJournalIds: Seq[Long]) = {
    dataAccess.summaryQueueEntries.filter(_.metadataJournalId.inSet(metadataJournalIds)).delete
  }

  private[slick] def countSummaryQueueEntries() = {
    dataAccess.summaryQueueEntries.length.result
  }

}
