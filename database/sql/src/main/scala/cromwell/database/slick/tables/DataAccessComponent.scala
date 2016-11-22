package cromwell.database.slick.tables

import slick.driver.JdbcProfile

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
    with CallCachingDetritusEntryComponent
    with CallCachingEntryComponent
    with CallCachingHashEntryComponent
    with CallCachingSimpletonEntryComponent
    with JobKeyValueEntryComponent
    with JobStoreEntryComponent
    with JobStoreSimpletonEntryComponent
    with MetadataEntryComponent
    with SummaryStatusEntryComponent
    with WorkflowMetadataSummaryEntryComponent
    with WorkflowStoreEntryComponent
    with SubWorkflowStoreEntryComponent {

  import driver.api._

  lazy val schema =
    callCachingDetritusEntries.schema ++
      callCachingEntries.schema ++
      callCachingHashEntries.schema ++
      callCachingSimpletonEntries.schema ++
      jobKeyValueEntries.schema ++
      jobStoreEntries.schema ++
      jobStoreSimpletonEntries.schema ++
      metadataEntries.schema ++
      summaryStatusEntries.schema ++
      workflowMetadataSummaryEntries.schema ++
      workflowStoreEntries.schema ++
      subWorkflowStoreEntries.schema
}
