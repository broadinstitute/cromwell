package cromwell.database.slick.tables

import slick.driver.JdbcProfile

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
    with MetadataComponent
    with WorkflowMetadataSummaryComponent
    with WorkflowStoreComponent
    with JobStoreComponent
    with JobStoreResultSimpletonComponent
    with BackendKVStoreComponent
    with CallCachingResultMetaInfoComponent
    with CallCachingHashComponent
    with CallCachingResultSimpletonComponent
    with SummaryStatusComponent
    with CallCachingJobDetritusComponent {

  import driver.api._

  lazy val schema =
    metadata.schema ++
      workflowMetadataSummaries.schema ++
      workflowStore.schema ++
      jobStore.schema ++
      jobStoreResultSimpletons.schema ++
      backendKVStore.schema ++
      callCachingResultMetaInfos.schema ++
      callCachingHashes.schema ++
      callCachingResultSimpletons.schema ++
      summaryStatuses.schema ++
      callCachingJobDetritus.schema
}
