package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingEntry

trait CallCachingEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class CallCachingEntries(tag: Tag) extends Table[CallCachingEntry](tag, "CALL_CACHING_ENTRY") {
    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME", O.Length(255))

    def jobIndex = column[Int]("JOB_INDEX")

    def jobAttempt = column[Option[Int]]("JOB_ATTEMPT", O.Default(None))

    def returnCode = column[Option[Int]]("RETURN_CODE")

    def allowResultReuse = column[Boolean]("ALLOW_RESULT_REUSE", O.Default(true))
    
    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, returnCode, allowResultReuse,
      callCachingEntryId.?) <> (CallCachingEntry.tupled, CallCachingEntry.unapply)

    def ucCallCachingEntryWeuCfqnJi =
      index("UC_CALL_CACHING_ENTRY_WEU_CFQN_JI", (workflowExecutionUuid, callFullyQualifiedName, jobIndex),
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

  val callCachingEntriesForWorkflowFqnIndex = Compiled(
    (workflowId: Rep[String], callFqn: Rep[String],  jobIndex: Rep[Int]) => for {
        callCachingEntry <- callCachingEntries
        if callCachingEntry.workflowExecutionUuid === workflowId
        if callCachingEntry.callFullyQualifiedName === callFqn
        if callCachingEntry.jobIndex === jobIndex
      } yield callCachingEntry
  )

  def callCachingEntryIdsForWorkflowId(workflowId: String) = {
    for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.workflowExecutionUuid === workflowId
    } yield callCachingEntry.callCachingEntryId
  }
}
