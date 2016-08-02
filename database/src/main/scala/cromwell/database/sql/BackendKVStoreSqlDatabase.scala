package cromwell.database.sql

import scala.concurrent.{ExecutionContext, Future}

trait BackendKVStoreSqlDatabase {
  this: SqlDatabase =>

  def addBackendStoreKeyValuePair(workflowUuid: String, callFqn: String, callIndex: Int, callAttempt: Int, storeKey: String, storeValue: String)(implicit ec: ExecutionContext): Future[Unit]

  def queryBackendStoreValueByStoreKey(workflowUid: String, callFqn: String, callIndex: Int, callAttempt: Int, jobKey: String)(implicit ec: ExecutionContext): Future[Option[String]]

}
