package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowExecution

@deprecated("Olde Worlde Databasee Tablee", "0.21")
trait WorkflowExecutionComponent {
  this: DriverComponent =>

  import driver.api._

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  class WorkflowExecutions(tag: Tag) extends Table[WorkflowExecution](tag, "WORKFLOW_EXECUTION") {
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def name = column[String]("WORKFLOW_NAME")

    def status = column[String]("STATUS")

    def startDt = column[Timestamp]("START_DT")

    def endDt = column[Option[Timestamp]]("END_DT")

    override def * = (workflowExecutionUuid, name, status, startDt, endDt, workflowExecutionId.?) <>
      (WorkflowExecution.tupled, WorkflowExecution.unapply)

    def uniqueKey = index("UK_WE_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, unique = true)

    def statusIdx = index("STATUS_IDX", status, unique = false)
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  protected val workflowExecutions = TableQuery[WorkflowExecutions]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionIdsAutoInc = workflowExecutions returning workflowExecutions.map(_.workflowExecutionId)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionIdsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid

    } yield workflowExecution.workflowExecutionId)
  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionStatusByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.status)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionStatusEndDtByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution.status, workflowExecution.endDt))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val workflowExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution)

}
