package cromwell.database.slick

import cromwell.database.sql.BackendKVStoreSqlDatabase
import cromwell.database.sql.tables.BackendKVStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait BackendKVStoreSlickDatabase extends BackendKVStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addBackendStoreKeyValuePair(workflowUuid: String, callFqn: String, callIndex: Int, callAttempt: Int, storeKey: String, storeValue: String)(implicit ec: ExecutionContext): Future[Unit] = {

    val entry = BackendKVStoreEntry(workflowUuid,
      callFqn,
      callIndex,
      callAttempt,
      storeKey,
      storeValue)

    val action = if (useSlickUpserts) {
      for {
        _ <- dataAccess.backendKVStoreAutoInc.insertOrUpdate(entry)
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          backendJobValueByBackendJobKey(workflowUuid, callFqn, callIndex, callAttempt, storeKey).
          update(storeValue)
        _ <- updateCount match {
          case 0 => dataAccess.backendKVStoreAutoInc += entry
          case _ => assertUpdateCount("addBackendJobKeyValuePair", updateCount, 1)
        }
      } yield () // <-- maps the result to a DBIO[Unit]
    }
    runTransaction(action)
  }

  override def queryBackendStoreValueByStoreKey(workflowUuid: String, callFqn: String, callIndex: Int, callAttempt: Int, jobKey: String)(implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.backendJobValueByBackendJobKey(workflowUuid, callFqn, callIndex, callAttempt, jobKey).result.headOption
    runTransaction(action)
  }
}
