package cromwell.database.slick.tables

import cromwell.database.sql.tables.JobStoreEntry

trait JobStoreComponent {

  this: DriverComponent =>

  import driver.api._

  class JobStoreEntries(tag: Tag) extends Table[JobStoreEntry](tag, "JOB_STORE") {
    def jobStoreId = column[Int]("JOB_STORE_ID", O.PrimaryKey, O.AutoInc)
    def workflowUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def callFqn = column[String]("CALL_FQN")
    def scatterIndex = column[Int]("JOB_SCATTER_INDEX")
    def attempt = column[Int]("JOB_RETRY_ATTEMPT")
    def jobSuccessful = column[Boolean]("JOB_SUCCESSFUL")
    def returnCode = column[Option[Int]]("RETURN_CODE")
    // Only set for failure:
    def exceptionMessage = column[Option[String]]("EXCEPTION_MESSAGE")
    def retryableFailure = column[Option[Boolean]]("RETRYABLE_FAILURE")


    override def * = (workflowUuid, callFqn, scatterIndex, attempt, jobSuccessful, returnCode, exceptionMessage, retryableFailure, jobStoreId.?) <>
      (JobStoreEntry.tupled, JobStoreEntry.unapply)

    def uuidIndex = index("JOB_STORE_UUID_IDX", workflowUuid, unique = false)
    def jobkeyIndex = index("JOB_STORE_JOBKEY_IDX", (workflowUuid, callFqn, scatterIndex, attempt), unique = true)
  }

  protected val jobStore = TableQuery[JobStoreEntries]

  val jobStoreAutoInc = jobStore returning jobStore.map(_.jobStoreId)

  /**
    * Useful for finding all store entry for a given workflow UUID (e.g. so you can delete them! Bwahaha)
    */
  val jobStoreEntryByWorkflowUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      jobStoreRow <- jobStore
      if jobStoreRow.workflowUuid === workflowExecutionUuid
    } yield jobStoreRow)

  /**
    * Useful for finding the unique job store entry for a given job key
    */
  val jobStoreEntriesByJobStoreKey = Compiled(
    (workflowUuid: Rep[String], callFqn: Rep[String], scatterIndex: Rep[Int], attempt: Rep[Int]) => for {
        jobStoreRow <- jobStore
        if jobStoreRow.workflowUuid === workflowUuid && jobStoreRow.callFqn === callFqn && jobStoreRow.scatterIndex === scatterIndex && jobStoreRow.attempt === attempt
      } yield jobStoreRow)
}
