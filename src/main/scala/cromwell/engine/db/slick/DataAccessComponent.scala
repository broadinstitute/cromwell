package cromwell.engine.db.slick

import slick.driver.JdbcProfile

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
  with WorkflowExecutionComponent
  with WorkflowExecutionAuxComponent
  with SymbolComponent
  with ExecutionComponent
  with JesJobComponent
  with LocalJobComponent
  with SgeJobComponent
  with ExecutionEventComponent {

  import driver.api._

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      localJobs.schema ++
      jesJobs.schema ++
      sgeJobs.schema ++
      executionEvents.schema
}
