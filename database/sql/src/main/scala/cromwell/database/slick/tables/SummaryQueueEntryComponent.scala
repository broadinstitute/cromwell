package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.SummaryQueueEntry

trait SummaryQueueEntryComponent {
  this: DriverComponent =>

  import driver.api._

  class SummaryQueueEntries(tag: Tag) extends Table[SummaryQueueEntry](tag, "SUMMARY_QUEUE_ENTRY") {

    def summaryQueueEntryId = column[Long]("SUMMARY_STATUS_ENTRY_ID", O.PrimaryKey, O.AutoInc)
    def metadataEntryId = column[Long]("METADATA_ENTRY_ID")
    def metadataEntryWriteTimestamp = column[Timestamp]("METADATA_ENTRY_WRITE_TIMESTAMP")
    def summarizationTimestamp = column[Option[Timestamp]]("SUMMARIZATION_TIMESTAMP")

    override def * = (metadataEntryId, metadataEntryWriteTimestamp, summarizationTimestamp, summaryQueueEntryId.?) <>
      (SummaryQueueEntry.tupled, SummaryQueueEntry.unapply)

    def ucsummaryQueueEntryMetadataEntryId = index("UC_SUMMARY_QUEUE_ENTRY_METADATA_ENTRY_ID", metadataEntryId, unique = true)
  }

  protected val summaryQueueEntries = TableQuery[SummaryQueueEntries]

  val summaryQueueEntryIdsAutoInc = summaryQueueEntries returning summaryQueueEntries.map(_.summaryQueueEntryId)

  val summaryQueueSummarizationTimestampForId = Compiled(
    (metadataEntryId: Rep[Long]) => for {
      entry <- summaryQueueEntries
      if entry.metadataEntryId === metadataEntryId
    } yield entry.summarizationTimestamp
  )

  def unsummarizedMetadataEntryIds(maxResults: Long) = (
    for {
      entry <- summaryQueueEntries
      if entry.summarizationTimestamp.isEmpty
    } yield entry.metadataEntryId).sortBy(identity).take(maxResults)

  def countUnsummarizedMetadataEntryIds() = {
    (for {
      entry <- summaryQueueEntries
      if entry.summarizationTimestamp.isEmpty
    } yield entry).length
  }
}

/**
val metadataEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )
  */
