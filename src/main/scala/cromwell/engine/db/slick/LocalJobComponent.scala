package cromwell.engine.db.slick

//TODO: Maybe replace the name of this to just, Job, say?
case class LocalJob
(
  executionId: Int,
  pid: Option[Int],
  localJobId: Option[Int] = None
  )

trait LocalJobComponent {
  this: DriverComponent with ExecutionComponent with WorkflowExecutionComponent =>

  import driver.api._

  class LocalJobs(tag: Tag) extends Table[LocalJob](tag, "LOCAL_JOB") {
    def localJobId = column[Int]("LOCAL_JOB_ID", O.PrimaryKey, O.AutoInc)
    def executionId = column[Int]("EXECUTION_ID")
    def pid = column[Option[Int]]("PID")

    override def * = (executionId, pid, localJobId.?) <>
      (LocalJob.tupled, LocalJob.unapply)

    def execution = foreignKey("FK_LOCAL_JOB_EXECUTION_ID", executionId, executions)(_.executionId)
    def uniqueKey = index("UK_LJ_EXECUTION_UUID", executionId, unique = true)
  }

  protected val localJobs = TableQuery[LocalJobs]

  val localJobsAutoInc = localJobs returning localJobs.
    map(_.localJobId) into ((a, id) => a.copy(localJobId = Some(id)))

  val localJobsByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      localJob <- localJobs
      if localJob.executionId === executionId
    } yield localJob)

  val localJobPidsByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      localJob <- localJobs
      if localJob.executionId === executionId
    } yield localJob.pid)

  val localJobsWithExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      execution <- executions
      if execution.workflowExecutionId === workflowExecution.workflowExecutionId
      localJob <- localJobs
      if localJob.executionId === execution.executionId
    } yield (execution, localJob))
}
