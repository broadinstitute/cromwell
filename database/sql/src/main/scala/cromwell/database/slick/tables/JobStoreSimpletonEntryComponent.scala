package cromwell.database.slick.tables

import cromwell.database.sql.tables.JobStoreSimpletonEntry
import slick.model.ForeignKeyAction.Cascade

trait JobStoreSimpletonEntryComponent {

  this: DriverComponent with JobStoreEntryComponent =>

  import driver.api._

  class JobStoreSimpletonEntries(tag: Tag) extends Table[JobStoreSimpletonEntry](tag, "JOB_STORE_SIMPLETON_ENTRY") {
    def jobStoreSimpletonEntryId = column[Int]("JOB_STORE_SIMPLETON_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def simpletonKey = column[String]("SIMPLETON_KEY")

    def simpletonValue = column[String]("SIMPLETON_VALUE")

    def wdlType = column[String]("WDL_TYPE")

    def jobStoreEntryId = column[Int]("JOB_STORE_ENTRY_ID")

    override def * = (simpletonKey, simpletonValue, wdlType, jobStoreEntryId.?, jobStoreSimpletonEntryId.?) <>
      (JobStoreSimpletonEntry.tupled, JobStoreSimpletonEntry.unapply)

    def fkJobStoreResultSimpletonJobStoreEntryId = foreignKey("FK_JOB_STORE_RESULT_SIMPLETON_JOB_STORE_ENTRY_ID",
      jobStoreEntryId, jobStoreEntries)(_.jobStoreEntryId, onDelete = Cascade)

    def ucJobStoreSimpletonEntryJseiSk =
      index("UC_JOB_STORE_SIMPLETON_ENTRY_JSEI_SK", (jobStoreEntryId, simpletonKey), unique = true)
  }

  protected val jobStoreSimpletonEntries = TableQuery[JobStoreSimpletonEntries]

  val jobStoreSimpletonEntryIdsAutoInc = jobStoreSimpletonEntries returning
    jobStoreSimpletonEntries.map(_.jobStoreSimpletonEntryId)

  /**
    * Find all result simpletons which match a given JOB_STORE_ENTRY_ID
    */
  val jobStoreSimpletonEntriesForJobStoreEntryId = Compiled(
    (jobStoreEntryId: Rep[Int]) => for {
      jobStoreSimpletonEntry <- jobStoreSimpletonEntries if jobStoreSimpletonEntry.jobStoreEntryId === jobStoreEntryId
    } yield jobStoreSimpletonEntry
  )
}
