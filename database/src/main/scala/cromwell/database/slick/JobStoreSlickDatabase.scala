package cromwell.database.slick

import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.tables.JobStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreSlickDatabase extends JobStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addJobStoreEntries(entries: Iterable[JobStoreEntry])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.jobStoreAutoInc ++= entries
    runTransaction(action) map { _ => () }
  }

  override def fetchJobResult(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                    (implicit ec: ExecutionContext): Future[Option[JobStoreEntry]] = {
    val action = dataAccess.jobStoreEntriesByJobStoreKey(workflowUuid, callFqn, index, attempt).result
    runTransaction(action) map {
      x => if (x.size == 1) Option(x.head)
      else if (x.size > 1) throw new Exception(s"That's impossible! Expected one JobResult for $workflowUuid:$callFqn:$index:$attempt but got ${x.size}")
      else None
    }
  }

  override def removeJobResultsForWorkflows(workflowUuids: Iterable[String])(implicit ec: ExecutionContext): Future[Int] = {
    val actions = workflowUuids map { workflowUuid =>  dataAccess.jobStoreEntryByWorkflowUuid(workflowUuid).delete }
    runTransaction(DBIO.sequence(actions)) map { _.sum }
  }
}
