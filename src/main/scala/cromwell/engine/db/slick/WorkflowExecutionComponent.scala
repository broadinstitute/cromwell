package cromwell.engine.db.slick

import java.sql.Timestamp
import java.util.UUID

import slick.dbio.{NoStream, Effect}

case class WorkflowExecution
(
  workflowExecutionUuid: String,
  status: String,
  startDt: Timestamp,
  workflowExecutionId: Option[Int] = None,
  endDt: Option[Timestamp] = None
  )

trait WorkflowExecutionComponent {
  this: DriverComponent =>

  import driver.api._

  class WorkflowExecutions(tag: Tag) extends Table[WorkflowExecution](tag, "WORKFLOW_EXECUTION") {
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def status = column[String]("STATUS")

    def startDt = column[Timestamp]("START_DT")

    def endDt = column[Option[Timestamp]]("END_DT")

    override def * = (workflowExecutionUuid, status, startDt, workflowExecutionId.?, endDt) <>
      (WorkflowExecution.tupled, WorkflowExecution.unapply)
  }

  val workflowExecutions = TableQuery[WorkflowExecutions]

  val workflowExecutionsAutoInc = workflowExecutions returning workflowExecutions.
    map(_.workflowExecutionId) into ((a, id) => a.copy(workflowExecutionId = Some(id)))

  val workflowExecutionByID = Compiled(
    (id: Rep[Int]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionId === id
    } yield workflowExecution)
}
