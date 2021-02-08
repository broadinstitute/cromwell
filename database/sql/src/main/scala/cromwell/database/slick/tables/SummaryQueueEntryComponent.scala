package cromwell.database.slick.tables

import cromwell.database.sql.tables.SummaryQueueEntry

trait SummaryQueueEntryComponent {
  this: DriverComponent =>

  import driver.api._

  class SummaryQueueEntries(tag: Tag) extends Table[SummaryQueueEntry](tag, "SUMMARY_QUEUE_ENTRY") {

    def metadataJournalId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey)

    override def * = metadataJournalId <> (SummaryQueueEntry.apply, SummaryQueueEntry.unapply)

  }

  val summaryQueueEntries = TableQuery[SummaryQueueEntries]

}
