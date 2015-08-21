package cromwell.engine.db.slick

import java.sql.Timestamp

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

    def uniqueKey = index("UK_WE_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, unique = true)
  }

  protected val workflowExecutions = TableQuery[WorkflowExecutions]

  val workflowExecutionsAutoInc = workflowExecutions returning workflowExecutions.
    map(_.workflowExecutionId) into ((a, id) => a.copy(workflowExecutionId = Some(id)))

  val workflowExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution)

  val workflowExecutionStatusesByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.status)

  // NOTE: No precompile for you!
  // [Compile] works for all functions ... consisting only of individual columns
  // - http://slick.typesafe.com/doc/3.0.0/queries.html
  // - https://groups.google.com/forum/#!topic/scalaquery/2d_r4DEthfY
  //
  // On one hand:
  // a) This should be compilable if we narrowed the API to only get Submitted/Running workflows... BUT
  // b) It turns out this is currently executed only once at server restart
  def workflowExecutionsByStatuses(statuses: Traversable[String]) = {
    for {
      workflowExecution <- workflowExecutions
      if workflowExecution.status inSet statuses
    } yield workflowExecution
  }
}
