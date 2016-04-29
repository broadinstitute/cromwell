package cromwell.database.slick

import slick.driver.JdbcProfile

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
  with WorkflowExecutionComponent
  with WorkflowExecutionAuxComponent
  with SymbolComponent
  with ExecutionComponent
  with ExecutionInfoComponent
  with ExecutionEventComponent
  with RuntimeAttributeComponent
  with FailureEventComponent {

  import driver.api._

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      executionInfos.schema ++
      executionEvents.schema ++
      runtimeAttributes.schema ++
      failureEvents.schema
}
