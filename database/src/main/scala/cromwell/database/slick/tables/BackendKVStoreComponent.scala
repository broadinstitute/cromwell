package cromwell.database.slick.tables

import cromwell.database.sql.tables.BackendKVStore

trait BackendKVStoreComponent {

  this: DriverComponent =>

  import driver.api._

  class BackendKeyValuePairs(tag: Tag) extends Table[BackendKVStore](tag, "BACKEND_KV_STORE") {
    def backendKVStoreID = column[Int]("BACKEND_KV_STORE_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def callFqn = column[String]("CALL_FQN")
    def callIndex = column[Int]("CALL_INDEX")
    def callAttempt = column[Int]("CALL_ATTEMPT")
    def backendJobKey = column[String]("BACKEND_JOB_KEY")
    def backendJobValue = column[String]("BACKEND_JOB_VALUE")


    override def * = (workflowExecutionUuid, callFqn, callIndex, callAttempt, backendJobKey, backendJobValue, backendKVStoreID.?) <>
      (BackendKVStore.tupled, BackendKVStore.unapply)

    def backendJobKeyIndex = index("UK_BACKEND_KV_STORE_JOB_KEY", (workflowExecutionUuid, callFqn, callIndex, callAttempt, backendJobKey), unique = true)
  }

  protected val backendKVStore = TableQuery[BackendKeyValuePairs]

  val backendKVStoreAutoInc = backendKVStore returning backendKVStore.map(_.backendKVStoreID)

  val backendJobValueByBackendJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], callIndex: Rep[Int], callAttempt: Rep[Int], backendJobKey: Rep[String]) => for {
      backendKeyValuePair <- backendKVStore
      if backendKeyValuePair.workflowExecutionUuid === workflowExecutionUuid
      if backendKeyValuePair.callFqn === callFqn
      if backendKeyValuePair.callIndex === callIndex
      if backendKeyValuePair.callAttempt === callAttempt
      if backendKeyValuePair.backendJobKey === backendJobKey
    } yield backendKeyValuePair.backendJobValue)
}
