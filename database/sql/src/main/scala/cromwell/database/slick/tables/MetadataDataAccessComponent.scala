package cromwell.database.slick.tables

import slick.jdbc.JdbcProfile

class MetadataDataAccessComponent(val driver: JdbcProfile) extends DataAccessComponent
  with CustomLabelEntryComponent
  with MetadataEntryComponent
  with SummaryStatusEntryComponent
  with WorkflowMetadataSummaryEntryComponent {

  import driver.api._

  override lazy val schema: driver.SchemaDescription =
      customLabelEntries.schema ++
      metadataEntries.schema ++
      summaryStatusEntries.schema ++
      workflowMetadataSummaryEntries.schema
}
