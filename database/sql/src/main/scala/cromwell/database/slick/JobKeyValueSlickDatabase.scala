package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.JobKeyValueSqlDatabase
import cromwell.database.sql.tables.JobKeyValueEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobKeyValueSlickDatabase extends JobKeyValueSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

  override def existsJobKeyValueEntries()(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.jobKeyValueEntriesExists.result
    runTransaction(action)
  }

  override def addJobKeyValueEntry(jobKeyValueEntry: JobKeyValueEntry)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = if (useSlickUpserts) {
      for {
        _ <- dataAccess.jobKeyValueEntryIdsAutoInc.insertOrUpdate(jobKeyValueEntry)
      } yield ()
    } else manualUpsertQuery(jobKeyValueEntry)
    runTransaction(action)
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!! Be careful using this function with multiple        !!!!!!!!
  // !!!!!!! updates running in a single transaction.            !!!!!!!!
  // !!!!!!! https://broadworkbench.atlassian.net/browse/BA-6262 !!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  private def manualUpsertQuery(jobKeyValueEntry: JobKeyValueEntry)
                       (implicit ec: ExecutionContext) = for {
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

  def addJobKeyValueEntries(jobKeyValueEntries: Iterable[JobKeyValueEntry])
                           (implicit ec: ExecutionContext): Future[Unit] = {
    val action = if (useSlickUpserts) {
      createBatchUpsert("KeyValueStore", dataAccess.jobKeyValueTableQueryCompiled, jobKeyValueEntries)
    } else {
      DBIO.sequence(jobKeyValueEntries.map(manualUpsertQuery))
    }
    runTransaction(action).void
  }

  override def queryJobKeyValueEntries(workflowExecutionUuid: String)
                                      (implicit ec: ExecutionContext): Future[Seq[JobKeyValueEntry]] = {
    val action = dataAccess.jobKeyValueEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result
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
