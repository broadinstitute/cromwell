package cromwell.database.slick.tables

import slick.driver.JdbcProfile

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
    with WorkflowExecutionComponent
    with WorkflowExecutionAuxComponent
    with SymbolComponent
    with ExecutionComponent
    with ExecutionInfoComponent
    with RuntimeAttributeComponent
    with MetadataComponent
    with WorkflowMetadataSummaryComponent
    with WorkflowStoreComponent
    with JobStoreComponent
    with JobStoreResultSimpletonComponent
    with BackendKVStoreComponent
    with CallCachingResultMetaInfoComponent
    with CallCachingHashComponent
    with CallCachingResultSimpletonComponent
    with SummaryStatusComponent {

  import driver.api._

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      executionInfos.schema ++
      runtimeAttributes.schema ++
      metadata.schema ++
      workflowMetadataSummaries.schema ++
      workflowStore.schema ++
      jobStore.schema ++
      jobStoreResultSimpletons.schema ++
      backendKVStore.schema ++
      callCachingResultMetaInfo.schema ++
      callCachingHashes.schema ++
      callCachingResultSimpletons.schema ++
      summaryStatuses.schema
}
