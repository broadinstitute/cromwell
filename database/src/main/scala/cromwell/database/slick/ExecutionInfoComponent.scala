package cromwell.database.slick

import cromwell.database.obj.ExecutionInfo

@deprecated("Olde Worlde Databasee Tablee", "0.21")
trait ExecutionInfoComponent {
  this: DriverComponent with ExecutionComponent with WorkflowExecutionComponent =>

  import driver.api._

  class ExecutionInfos(tag: Tag) extends Table[ExecutionInfo](tag, "EXECUTION_INFO") {
    def executionInfoId = column[Int]("EXECUTION_INFO_ID", O.PrimaryKey, O.AutoInc)
    def executionId = column[Int]("EXECUTION_ID")
    def key = column[String]("INFO_KEY")
    def value = column[Option[String]]("INFO_VALUE")
    override def * = (executionId, key, value, executionInfoId.?) <> (ExecutionInfo.tupled, ExecutionInfo.unapply)
    def execution = foreignKey("FK_EXECUTION_INFO_EXECUTION", executionId, executions)(_.executionId)
    def uniqueKey = index("UK_EXECUTION_INFO", (executionId, key), unique = true)
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  protected val executionInfos = TableQuery[ExecutionInfos]
  val executionInfoIdsAutoInc = executionInfos returning executionInfos.map(_.executionInfoId)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionInfosByWorkflowExecutionUuidAndCallFqnAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], attempt: Rep[Int]) => for {
      executionInfo <- executionInfos
      execution <- executionInfo.execution
      if execution.callFqn === callFqn
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield executionInfo)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionInfoValueByExecutionAndKey = Compiled(
    (executionId: Rep[Int], key: Rep[String]) => for {
      executionInfo <- executionInfos
      if executionInfo.key === key
      if executionInfo.executionId === executionId
    } yield executionInfo.value)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionInfoValueByWorkflowExecutionUuidAndCallFqnAndAttemptAndKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], attempt: Rep[Int], key: Rep[String]) => for {
      executionInfo <- executionInfos
      if executionInfo.key === key
      execution <- executionInfo.execution
      if execution.callFqn === callFqn
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield executionInfo.value)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionsAndExecutionInfosByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      executionInfo <- executionInfos
      execution <- executionInfo.execution
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution, executionInfo)
  )

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionsAndExecutionInfosByWorkflowExecutionUuidAndCallFqn = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String]) => for {
      executionInfo <- executionInfos
      execution <- executionInfo.execution
      if execution.callFqn === callFqn
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution, executionInfo)
  )

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def runningExecutionsAndExecutionInfosByWorkflowExecutionUuid(workflowExecutionUuid: String,
                                                                statuses: Traversable[String]) = {
    for {
      executionInfo <- executionInfos
      execution <- executionInfo.execution
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      if execution.status inSet statuses
    } yield (execution, executionInfo)
  }
}
