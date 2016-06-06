package cromwell.database.slick

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
  with WorkflowMetadataSummaryComponent {

  import driver.api._

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      executionInfos.schema ++
      runtimeAttributes.schema ++
      metadata.schema ++
      workflowMetadataSummaries.schema
}
