package cromwell.engine.db.slick

case class ExecutionInfo(executionId: Int,
                         key: String,
                         value: Option[String] = None,
                         executionInfoId: Option[Int] = None)

case class ExecutionAndExecutionInfo(execution: Execution, executionInfo: ExecutionInfo)
case class ExecutionInfosByExecution(execution: Execution, executionInfos: Seq[ExecutionInfo])

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

  protected val executionInfos = TableQuery[ExecutionInfos]
  val executionInfosAutoInc = executionInfos returning executionInfos.map(_.executionInfoId) into ((a, id) => a.copy(executionInfoId = Option(id)))

  val executionInfosByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      executionInfo <- executionInfos
      if executionInfo.executionId === executionId
    } yield executionInfo)

  val executionInfoValueByExecutionAndKey = Compiled(
    (executionId: Rep[Int], key: Rep[String]) => for {
      executionInfo <- executionInfos
      if executionInfo.key === key
      if executionInfo.executionId === executionId
    } yield executionInfo.value)

  val executionsAndExecutionInfosByWorkflowId = Compiled(
    (id: Rep[String]) => for {
      executionInfo <- executionInfos
      execution <- executionInfo.execution
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === id
    } yield (execution, executionInfo)
  )
}
