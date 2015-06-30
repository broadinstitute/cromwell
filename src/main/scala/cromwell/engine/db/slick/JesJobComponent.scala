package cromwell.engine.db.slick

case class JesJob
(
  executionId: Int,
  jesId: Int,
  jesStatus: String,
  jesJobId: Option[Int] = None
  )

trait JesJobComponent {
  this: DriverComponent with ExecutionComponent =>

  import driver.api._

  class JesJobs(tag: Tag) extends Table[JesJob](tag, "JES_JOB") {
    def jesJobId = column[Int]("JES_JOB_ID", O.PrimaryKey, O.AutoInc)

    def executionId = column[Int]("EXECUTION_ID")

    def jesId = column[Int]("JES_ID")

    def jesStatus = column[String]("JES_STATUS")

    override def * = (executionId, jesId, jesStatus, jesJobId.?) <>
      (JesJob.tupled, JesJob.unapply)

    def execution = foreignKey(
      "FK_JES_JOB_EXECUTION_ID", executionId, executions)(_.executionId)

    def uniqueKey = index("UK_JJ_EXECUTION_UUID",
      executionId, unique = true)
  }

  val jesJobs = TableQuery[JesJobs]

  val jesJobsAutoInc = jesJobs returning jesJobs.
    map(_.jesJobId) into ((a, id) => a.copy(jesJobId = Some(id)))

  val jesJobByID = Compiled(
    (id: Rep[Int]) => for {
      jesJob <- jesJobs
      if jesJob.jesJobId === id
    } yield jesJob)
}
