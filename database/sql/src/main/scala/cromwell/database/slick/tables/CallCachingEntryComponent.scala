package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingEntry
import slick.profile.RelationalProfile.ColumnOption.Default

trait CallCachingEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class CallCachingEntries(tag: Tag) extends Table[CallCachingEntry](tag, "CALL_CACHING_ENTRY") {
    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME")

    def jobIndex = column[Int]("JOB_INDEX")

    def returnCode = column[Option[Int]]("RETURN_CODE")

    def allowResultReuse = column[Boolean]("ALLOW_RESULT_REUSE", Default(true))

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, returnCode, allowResultReuse,
      callCachingEntryId.?) <> (CallCachingEntry.tupled, CallCachingEntry.unapply)

    def ucCallCachingEntryWeuCqfnJi =
      index("UC_CALL_CACHING_ENTRY_WEU_CQFN_JI", (workflowExecutionUuid, callFullyQualifiedName, jobIndex),
        unique = true)
  }

  protected val callCachingEntries = TableQuery[CallCachingEntries]

  val callCachingEntryIdsAutoInc = callCachingEntries returning callCachingEntries.map(_.callCachingEntryId)

  val callCachingEntriesForId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingEntry
  )

  val allowResultReuseForCallCachingEntryId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingEntry.allowResultReuse
  )
}
