package cromwell.database.sql

import cromwell.database.sql.tables.BackendKVStore
import scala.concurrent.{ExecutionContext, Future}

trait BackendKVStoreSqlDatabase {
  this: SqlDatabase =>

  def addBackendJobKeyValuePair(entries: Iterable[BackendKVStore])(implicit ec: ExecutionContext): Future[Unit]

  def queryBackendJobValueByJobKey(workflowUid: String, callFqn: String, callIndex: Int, callAttempt: Int, jobKey: String)(implicit ec: ExecutionContext): Future[Option[String]]

}