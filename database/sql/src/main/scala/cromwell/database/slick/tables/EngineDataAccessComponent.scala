package cromwell.database.slick.tables

import slick.jdbc.JdbcProfile

class EngineDataAccessComponent(val driver: JdbcProfile)
    extends DataAccessComponent
    with CallCachingDetritusEntryComponent
    with CallCachingEntryComponent
    with CallCachingHashEntryComponent
    with CallCachingAggregationEntryComponent
    with CallCachingSimpletonEntryComponent
    with DockerHashStoreEntryComponent
    with JobKeyValueEntryComponent
    with JobStoreEntryComponent
    with JobStoreSimpletonEntryComponent
    with SubWorkflowStoreEntryComponent
    with WorkflowStoreEntryComponent {

  import driver.api._

  override lazy val schema: driver.SchemaDescription =
    callCachingDetritusEntries.schema ++
      callCachingEntries.schema ++
      callCachingHashEntries.schema ++
      callCachingSimpletonEntries.schema ++
      callCachingAggregationEntries.schema ++
      dockerHashStoreEntries.schema ++
      jobKeyValueEntries.schema ++
      jobStoreEntries.schema ++
      jobStoreSimpletonEntries.schema ++
      subWorkflowStoreEntries.schema ++
      workflowStoreEntries.schema
}
