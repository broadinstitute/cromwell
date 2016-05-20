package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.obj.WorkflowExecution

trait WorkflowExecutionComponent {
  this: DriverComponent =>

  import driver.api._

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

  protected val workflowExecutions = TableQuery[WorkflowExecutions]

  val workflowExecutionIdsAutoInc = workflowExecutions returning workflowExecutions.map(_.workflowExecutionId)

  val workflowExecutionIdsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.workflowExecutionId)

  val workflowExecutionStatusByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.status)

  val workflowExecutionStatusEndDtByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution.status, workflowExecution.endDt))

  val workflowExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution)

}
