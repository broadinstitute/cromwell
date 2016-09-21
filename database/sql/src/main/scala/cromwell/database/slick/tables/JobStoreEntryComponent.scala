package cromwell.database.slick.tables

import cromwell.database.sql.tables.JobStoreEntry

trait JobStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class JobStoreEntries(tag: Tag) extends Table[JobStoreEntry](tag, "JOB_STORE_ENTRY") {
    def jobStoreEntryId = column[Int]("JOB_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME")

    def jobIndex = column[Int]("JOB_INDEX")

    def jobAttempt = column[Int]("JOB_ATTEMPT")

    def jobSuccessful = column[Boolean]("JOB_SUCCESSFUL")

    def returnCode = column[Option[Int]]("RETURN_CODE")

    // Only set for failure:
    def exceptionMessage = column[Option[String]]("EXCEPTION_MESSAGE")

    def retryableFailure = column[Option[Boolean]]("RETRYABLE_FAILURE")

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, jobSuccessful, returnCode,
      exceptionMessage, retryableFailure, jobStoreEntryId.?) <> (JobStoreEntry.tupled, JobStoreEntry.unapply)

    def ucJobStoreEntryWeuCfqnJiJa = index("UC_JOB_STORE_ENTRY_WEU_CFQN_JI_JA",
      (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt), unique = true)

    def ixJobStoreEntryWeu = index("IX_JOB_STORE_ENTRY_WEU", workflowExecutionUuid, unique = false)
  }

  protected val jobStoreEntries = TableQuery[JobStoreEntries]

  val jobStoreEntryIdsAutoInc = jobStoreEntries returning jobStoreEntries.map(_.jobStoreEntryId)

  /**
    * Useful for finding all job stores for a given workflow execution UUID (e.g. so you can delete them! Bwahaha)
    */
  val jobStoreEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      jobStoreEntry <- jobStoreEntries
      if jobStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield jobStoreEntry
  )

  /**
    * Useful for finding the unique job store for a given job key
    */
  val jobStoreEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Int],
     jobAttempt: Rep[Int]) =>
      for {
        jobStoreEntry <- jobStoreEntries
        if jobStoreEntry.workflowExecutionUuid === workflowExecutionUuid &&
          jobStoreEntry.callFullyQualifiedName === callFullyQualifiedName &&
          jobStoreEntry.jobIndex === jobIndex && jobStoreEntry.jobAttempt === jobAttempt
      } yield jobStoreEntry
  )
}
