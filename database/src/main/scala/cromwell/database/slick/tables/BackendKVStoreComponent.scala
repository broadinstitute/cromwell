package cromwell.database.slick.tables

import cromwell.database.sql.tables.BackendKVStoreEntry

trait BackendKVStoreComponent {

  this: DriverComponent =>

  import driver.api._

  class BackendKeyValuePairs(tag: Tag) extends Table[BackendKVStoreEntry](tag, "BACKEND_KV_STORE") {
    def backendKVStoreID = column[Int]("BACKEND_KV_STORE_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def callFqn = column[String]("CALL_FQN")
    def jobScatterIndex = column[Int]("JOB_SCATTER_INDEX")
    def jobRetryAttempt = column[Int]("JOB_RETRY_ATTEMPT")
    def storeKey = column[String]("STORE_KEY")
    def storeValue = column[String]("STORE_VALUE")

    override def * = (workflowExecutionUuid, callFqn, jobScatterIndex, jobRetryAttempt, storeKey, storeValue, backendKVStoreID.?) <>
      (BackendKVStoreEntry.tupled, BackendKVStoreEntry.unapply)

    def backendKVStoreIndex = index("UK_BACKEND_KV_STORE_KEY", (workflowExecutionUuid, callFqn, jobScatterIndex, jobRetryAttempt, storeKey), unique = true)
  }

  protected val backendKVStore = TableQuery[BackendKeyValuePairs]

  val backendKVStoreAutoInc = backendKVStore returning backendKVStore.map(_.backendKVStoreID)

  val backendJobValueByBackendJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], jobScatterIndex: Rep[Int], jobRetryAttempt: Rep[Int], storeKey: Rep[String]) => for {
      backendKeyValuePair <- backendKVStore
      if backendKeyValuePair.workflowExecutionUuid === workflowExecutionUuid
      if backendKeyValuePair.callFqn === callFqn
      if backendKeyValuePair.jobScatterIndex === jobScatterIndex
      if backendKeyValuePair.jobRetryAttempt === jobRetryAttempt
      if backendKeyValuePair.storeKey === storeKey
    } yield backendKeyValuePair.storeValue)
}
