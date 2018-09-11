package cromwell.database.slick.tables

import cromwell.database.sql.tables.JobKeyValueEntry

trait JobKeyValueEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class JobKeyValueEntries(tag: Tag) extends Table[JobKeyValueEntry](tag, "JOB_KEY_VALUE_ENTRY") {
    def jobKeyValueEntryId = column[Int]("JOB_KEY_VALUE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME", O.Length(255))

    def jobIndex = column[Int]("JOB_INDEX")

    def jobAttempt = column[Int]("JOB_ATTEMPT")

    def storeKey = column[String]("STORE_KEY", O.Length(255))

    def storeValue = column[String]("STORE_VALUE", O.Length(255))

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, storeKey, storeValue,
      jobKeyValueEntryId.?) <> (JobKeyValueEntry.tupled, JobKeyValueEntry.unapply)

    def ucJobKeyValueEntryWeuCfqnJiJaSk = index("UC_JOB_KEY_VALUE_ENTRY_WEU_CFQN_JI_JA_SK",
      (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, storeKey), unique = true)
  }

  protected val jobKeyValueEntries = TableQuery[JobKeyValueEntries]
  lazy val jobKeyValueTableQueryCompiled = driver.compileInsert(jobKeyValueEntries.toNode)

  val jobKeyValueEntryIdsAutoInc = jobKeyValueEntries returning jobKeyValueEntries.map(_.jobKeyValueEntryId)

  val jobKeyValueEntriesForWorkflowExecutionUuid = Compiled((workflowExecutionUuid: Rep[String]) => for {
      jobKeyValueEntry <- jobKeyValueEntries
      if jobKeyValueEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield jobKeyValueEntry
  )

  val storeValuesForJobKeyAndStoreKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Int], jobAttempt: Rep[Int],
     storeKey: Rep[String]) => for {
      jobKeyValueEntry <- jobKeyValueEntries
      if jobKeyValueEntry.workflowExecutionUuid === workflowExecutionUuid
      if jobKeyValueEntry.callFullyQualifiedName === callFullyQualifiedName
      if jobKeyValueEntry.jobIndex === jobIndex
      if jobKeyValueEntry.jobAttempt === jobAttempt
      if jobKeyValueEntry.storeKey === storeKey
    } yield jobKeyValueEntry.storeValue
  )
}
