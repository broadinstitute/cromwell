package cromwell.engine.db.slick

case class LocalJob
(
  executionId: Int,
  pid: Int,
  rc: Int,
  localJobId: Option[Int] = None
  )

trait LocalJobComponent {
  this: DriverComponent with ExecutionComponent =>

  import driver.api._

  class LocalJobs(tag: Tag) extends Table[LocalJob](tag, "LOCAL_JOB") {
    def localJobId = column[Int]("LOCAL_JOB_ID", O.PrimaryKey, O.AutoInc)

    def executionId = column[Int]("EXECUTION_ID")

    def pid = column[Int]("PID")

    def rc = column[Int]("RC")

    override def * = (executionId, pid, rc, localJobId.?) <>
      (LocalJob.tupled, LocalJob.unapply)

    def execution = foreignKey(
      "FK_LOCAL_JOB_EXECUTION_ID", executionId, executions)(_.executionId)
  }

  val localJobs = TableQuery[LocalJobs]

  val localJobsAutoInc = localJobs returning localJobs.
    map(_.localJobId) into ((a, id) => a.copy(localJobId = Some(id)))

  val localJobByID = Compiled(
    (id: Rep[Int]) => for {
      localJob <- localJobs
      if localJob.localJobId === id
    } yield localJob)
}
