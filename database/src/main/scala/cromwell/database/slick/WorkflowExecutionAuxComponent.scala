package cromwell.database.slick

import java.sql.Clob

import cromwell.database.obj.WorkflowExecutionAux

trait WorkflowExecutionAuxComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

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

  protected val workflowExecutionAuxes = TableQuery[WorkflowExecutionAuxes]

  val workflowExecutionAuxIdsAutoInc = workflowExecutionAuxes returning workflowExecutionAuxes.
    map(_.workflowExecutionAuxId)

  val workflowExecutionsAndAuxesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.workflowExecutionId === workflowExecutionId
    } yield (workflowExecution, workflowExecutionAux))

  val workflowExecutionsAndAuxesByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution, workflowExecutionAux))

  val workflowOptionsByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      workflowExecutionAux <- workflowExecutionAuxes
      if workflowExecutionAux.workflowExecutionId === workflowExecutionId
    } yield workflowExecutionAux.workflowOptions)

  // NOTE: No precompile for you!
  // [Compile] works for all functions ... consisting only of individual columns
  // - http://slick.typesafe.com/doc/3.0.0/queries.html
  // - https://groups.google.com/forum/#!topic/scalaquery/2d_r4DEthfY
  def workflowExecutionAndAuxesByStatuses(statuses: Traversable[String]) = {
    for {
      workflowExecutionAux <- workflowExecutionAuxes
      workflowExecution <- workflowExecutionAux.workflowExecution
      if workflowExecution.status inSet statuses
    } yield (workflowExecution, workflowExecutionAux)
  }
}
