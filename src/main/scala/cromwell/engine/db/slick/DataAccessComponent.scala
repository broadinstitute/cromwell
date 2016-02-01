package cromwell.engine.db.slick

import slick.driver.JdbcProfile

//TODO: Remove backend references from here
class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
  with WorkflowExecutionComponent
  with WorkflowExecutionAuxComponent
  with SymbolComponent
  with ExecutionComponent
  with LocalJobComponent
  with ExecutionEventComponent {

  import driver.api._

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      localJobs.schema ++
      executionEvents.schema
}
