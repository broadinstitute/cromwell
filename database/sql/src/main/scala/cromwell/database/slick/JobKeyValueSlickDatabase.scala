package cromwell.database.slick

import cromwell.database.sql.JobKeyValueSqlDatabase
import cromwell.database.sql.tables.JobKeyValueEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobKeyValueSlickDatabase extends JobKeyValueSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addJobKeyValueEntry(jobKeyValueEntry: JobKeyValueEntry)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = if (useSlickUpserts) {
      for {
        _ <- dataAccess.jobKeyValueEntryIdsAutoInc.insertOrUpdate(jobKeyValueEntry)
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          storeValuesForJobKeyAndStoreKey((
            jobKeyValueEntry.workflowExecutionUuid,
            jobKeyValueEntry.callFullyQualifiedName,
            jobKeyValueEntry.jobIndex,
            jobKeyValueEntry.jobAttempt,
            jobKeyValueEntry.storeKey)).
          update(jobKeyValueEntry.storeValue)
        _ <- updateCount match {
          case 0 => dataAccess.jobKeyValueEntryIdsAutoInc += jobKeyValueEntry
          case _ => assertUpdateCount("addJobKeyValueEntry", updateCount, 1)
        }
      } yield ()
    }
    runTransaction(action)
  }

  override def queryStoreValue(workflowExecutionUuid: String, callFqn: String, jobScatterIndex: Int,
                               jobRetryAttempt: Int, storeKey: String)
                              (implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.
      storeValuesForJobKeyAndStoreKey((workflowExecutionUuid, callFqn, jobScatterIndex, jobRetryAttempt, storeKey)).
      result.headOption
    runTransaction(action)
  }
}
