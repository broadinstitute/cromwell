package cromwell.engine.db.slick

case class JesJob
(
  executionId: Int,
  jesId: Option[String],
  jesStatus: Option[String],
  jesJobId: Option[Int] = None
  )

trait JesJobComponent {
  this: DriverComponent with ExecutionComponent with WorkflowExecutionComponent =>

  import driver.api._

  class JesJobs(tag: Tag) extends Table[JesJob](tag, "JES_JOB") {
    def jesJobId = column[Int]("JES_JOB_ID", O.PrimaryKey, O.AutoInc)
    def executionId = column[Int]("EXECUTION_ID")
    def jesId = column[Option[String]]("JES_ID")
    def jesStatus = column[Option[String]]("JES_STATUS")

    override def * = (executionId, jesId, jesStatus, jesJobId.?) <>
      (JesJob.tupled, JesJob.unapply)

    def execution = foreignKey(
      "FK_JES_JOB_EXECUTION_ID", executionId, executions)(_.executionId)

    def uniqueKey = index("UK_JJ_EXECUTION_UUID",
      executionId, unique = true)
  }

  protected val jesJobs = TableQuery[JesJobs]

  val jesJobsAutoInc = jesJobs returning jesJobs.
    map(_.jesJobId) into ((a, id) => a.copy(jesJobId = Some(id)))

  val jesJobsByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      jesJob <- jesJobs
      if jesJob.executionId === executionId
    } yield jesJob)

  val jesIdsAndJesStatusesByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      jesJob <- jesJobs
      if jesJob.executionId === executionId
    } yield (jesJob.jesId, jesJob.jesStatus))

  val jesJobsWithExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
    workflowExecution <- workflowExecutions if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    execution <- executions if execution.workflowExecutionId === workflowExecution.workflowExecutionId
    jesJob <- jesJobs if jesJob.executionId === execution.executionId
  } yield (execution, jesJob))
}
