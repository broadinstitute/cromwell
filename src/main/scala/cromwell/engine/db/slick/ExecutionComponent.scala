package cromwell.engine.db.slick

case class Execution
(
  workflowExecutionId: Int,
  callFqn: String,
  status: String,
  executionId: Option[Int] = None
  )

trait ExecutionComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  class Executions(tag: Tag) extends Table[Execution](tag, "EXECUTION") {
    def executionId = column[Int]("EXECUTION_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")

    def callFqn = column[String]("CALL_FQN")

    def status = column[String]("STATUS")

    override def * = (workflowExecutionId, callFqn, status, executionId.?) <>
      (Execution.tupled, Execution.unapply)

    def workflow = foreignKey(
      "FK_EXECUTION_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)
  }

  val executions = TableQuery[Executions]

  val executionsAutoInc = executions returning executions.
    map(_.executionId) into ((a, id) => a.copy(executionId = Some(id)))

  val executionByID = Compiled(
    (id: Rep[Int]) => for {
      execution <- executions
      if execution.executionId === id
    } yield execution)
}
