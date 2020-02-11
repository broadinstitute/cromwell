package cromwell.database.slick.tables

import slick.jdbc.JdbcProfile

class MetadataDataAccessComponent(val driver: JdbcProfile) extends DataAccessComponent
  with CustomLabelEntryComponent
  with MetadataEntryComponent
  with SummaryStatusEntryComponent
  with SummaryQueueEntryComponent
  with WorkflowMetadataSummaryEntryComponent {

  import driver.api._

  override lazy val schema: driver.SchemaDescription =
      customLabelEntries.schema ++
      metadataEntries.schema ++
      summaryStatusEntries.schema ++
      workflowMetadataSummaryEntries.schema ++
      summaryQueueEntries.schema

  // Looks like here is the most appropriate place for this val since it doesn't fit neither in
  // SummaryQueueEntryComponent nor in MetadataEntryComponent
  val metadataEntriesToSummarizeQuery = {
    Compiled(
      (limit: ConstColumn[Long]) => (for {
        (_, metadataEntry) <-
          summaryQueueEntries join metadataEntries on (_.metadataJournalId === _.metadataEntryId)
      } yield metadataEntry).sortBy(_.metadataEntryId).take(limit)
    )
  }
}
