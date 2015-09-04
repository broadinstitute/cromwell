package cromwell.engine.db.slick

case class SgeJob (executionId: Int,
                   sgeJobNumber: Int,
                   sgeJobId: Option[Int] = None)

trait SgeJobComponent {
  this: DriverComponent with ExecutionComponent =>

  import driver.api._

  class SgeJobs(tag: Tag) extends Table[SgeJob](tag, "SGE_JOB") {
    def sgeJobId = column[Int]("SGE_JOB_ID", O.PrimaryKey, O.AutoInc)

    def executionId = column[Int]("EXECUTION_ID")

    def sgeJobNumber = column[Int]("SGE_JOB_NUMBER")

    override def * = (executionId, sgeJobNumber, sgeJobId.?) <> (SgeJob.tupled, SgeJob.unapply)

    def execution = foreignKey("FK_SGE_JOB_EXECUTION_ID", executionId, executions)(_.executionId)

    def uniqueKey = index("UK_SGE_JOB_EXECUTION_UUID", executionId, unique = true)
  }

  protected val sgeJobs = TableQuery[SgeJobs]

  val sgeJobsAutoInc = sgeJobs returning sgeJobs.
    map(_.sgeJobId) into ((a, id) => a.copy(sgeJobId = Some(id)))

  val sgeJobsByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      sgeJob <- sgeJobs
      if sgeJob.executionId === executionId
    } yield sgeJob)

  val sgeJobNumberByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      sgeJob <- sgeJobs
      if sgeJob.executionId === executionId
    } yield sgeJob.sgeJobNumber)

}
