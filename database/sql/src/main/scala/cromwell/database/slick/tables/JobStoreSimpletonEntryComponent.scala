package cromwell.database.slick.tables

import cromwell.database.sql.tables.JobStoreResultSimpletonEntry
import slick.model.ForeignKeyAction.Cascade


trait JobStoreResultSimpletonComponent {

  this: DriverComponent with JobStoreComponent =>

  import driver.api._

  class JobStoreResultSimpletonEntries(tag: Tag) extends Table[JobStoreResultSimpletonEntry](tag, "JOB_STORE_RESULT_SIMPLETON") {
    def jobStoreSimpletonId = column[Int]("JOB_STORE_RESULT_SIMPLETON_ID", O.PrimaryKey, O.AutoInc)
    def simpletonKey = column[String]("SIMPLETON_KEY")
    def simpletonValue = column[String]("SIMPLETON_VALUE")
    def wdlType = column[String]("WDL_TYPE")
    def jobStoreId = column[Int]("JOB_STORE_ID")

    override def * = (simpletonKey, simpletonValue, wdlType, jobStoreId, jobStoreSimpletonId.?) <>
      (JobStoreResultSimpletonEntry.tupled, JobStoreResultSimpletonEntry.unapply)

    def jobStoreResultSimpletonsUniquenessConstraint = index("UK_JOB_STORE_RESULT_SIMPLETON", (simpletonKey, jobStoreId), unique = true)
    def jobStoreForeignKey = foreignKey("JSRS_JOB_STORE_FK", jobStoreId, jobStore)(_.jobStoreId, onDelete = Cascade)
  }

  protected val jobStoreResultSimpletons = TableQuery[JobStoreResultSimpletonEntries]

  val jobStoreResultSimpletonAutoInc = jobStoreResultSimpletons returning jobStoreResultSimpletons.map(_.jobStoreSimpletonId)

  /**
    * Find all result simpletons which match a given JOB_STORE_ID
    */
  val jobStoreResultSimpletonsForJobStoreId = Compiled(
    (jobStoreId: Rep[Int]) => for {
      simpleton <- jobStoreResultSimpletons if simpleton.jobStoreId === jobStoreId
    } yield simpleton)
}

