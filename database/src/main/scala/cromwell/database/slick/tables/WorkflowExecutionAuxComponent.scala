package cromwell.database.slick.tables

import java.sql.Clob

import cromwell.database.sql.tables.WorkflowExecutionAux

@deprecated("Olde Worlde Databasee Tablee", "0.21")
trait WorkflowExecutionAuxComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  class WorkflowExecutionAuxes(tag: Tag) extends Table[WorkflowExecutionAux](tag, "WORKFLOW_EXECUTION_AUX") {
    def workflowExecutionAuxId = column[Int]("WORKFLOW_EXECUTION_AUX_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")
    def wdlSource = column[Clob]("WDL_SOURCE")
    def jsonInputs = column[Clob]("JSON_INPUTS")
    def workflowOptions = column[Clob]("WORKFLOW_OPTIONS")

    override def * = (workflowExecutionId, wdlSource, jsonInputs, workflowOptions, workflowExecutionAuxId.?) <>
      (WorkflowExecutionAux.tupled, WorkflowExecutionAux.unapply)

    def workflowExecution = foreignKey(
      "FK_WE_AUX_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_WE_AUX_WORKFLOW_EXECUTION_ID",
      workflowExecutionId, unique = true)
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  protected val workflowExecutionAuxes = TableQuery[WorkflowExecutionAuxes]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionAuxIdsAutoInc = workflowExecutionAuxes returning workflowExecutionAuxes.
    map(_.workflowExecutionAuxId)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionsAndAuxesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.workflowExecutionId === workflowExecutionId
    } yield (workflowExecution, workflowExecutionAux))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionsAndAuxesByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution, workflowExecutionAux))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowOptionsByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      if workflowExecutionAux.workflowExecutionId === workflowExecutionId
    } yield workflowExecutionAux.workflowOptions)

  // NOTE: No precompile for you!
  // [Compile] works for all functions ... consisting only of individual columns
  // - http://slick.typesafe.com/doc/3.0.0/queries.html
  // - https://groups.google.com/forum/#!topic/scalaquery/2d_r4DEthfY
  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def workflowExecutionAndAuxesByStatuses(statuses: Traversable[String]) = {
    for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.status inSet statuses
    } yield (workflowExecution, workflowExecutionAux)
  }
}
