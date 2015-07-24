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

    def workflowExecution = foreignKey(
      "FK_EXECUTION_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_EX_WORKFLOW_EXECUTION_ID",
      (workflowExecutionId, callFqn), unique = true)
  }

  protected val executions = TableQuery[Executions]

  val executionsAutoInc = executions returning executions.
    map(_.executionId) into ((a, id) => a.copy(executionId = Some(id)))

  val executionCallFqnsAndStatusesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield (execution.callFqn, execution.status))

  val executionStatusByWorkflowExecutionIdAndCallFqn = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
    } yield execution.status)

  val executionsByWorkflowExecutionUuidAndCallFqn = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  val executionStatusesByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      execution <- executions
      if execution.executionId === executionId
    } yield execution.status)

  // see workflowExecutionsByStatuses
  def executionStatusesByWorkflowExecutionIdAndCallFqns(workflowExecutionId: Int, callFqns: Traversable[String]) = {
    for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId && (execution.callFqn inSet callFqns)
    } yield execution.status
  }
}
