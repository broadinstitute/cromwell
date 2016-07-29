package cromwell.database.slick

import cromwell.database.sql.BackendKVStoreSqlDatabase
import cromwell.database.sql.tables.BackendKVStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait BackendKVStoreSlickDatabase extends BackendKVStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._


  override def addBackendKVStoreEntry(entry: BackendKVStoreEntry)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.backendKVStoreAutoInc += entry
    runTransaction(action) map { _ => () }
  }

  override def queryBackendKVStoreValueByStoreKey(workflowUid: String, callFqn: String, callIndex: Int, callAttempt: Int, jobKey: String)(implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.backendJobValueByBackendJobKey(workflowUid, callFqn, callIndex, callAttempt, jobKey).result.headOption
    runTransaction(action)
  }
}

