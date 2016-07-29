package cromwell.database.sql

import cromwell.database.sql.tables.BackendKVStoreEntry
import scala.concurrent.{ExecutionContext, Future}

trait BackendKVStoreSqlDatabase {
  this: SqlDatabase =>

  def addBackendKVStoreEntry(entry: BackendKVStoreEntry)(implicit ec: ExecutionContext): Future[Unit]

  def queryBackendKVStoreValueByStoreKey(workflowUid: String, callFqn: String, callIndex: Int, callAttempt: Int, jobKey: String)(implicit ec: ExecutionContext): Future[Option[String]]

}
