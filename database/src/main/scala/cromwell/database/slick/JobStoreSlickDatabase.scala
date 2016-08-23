package cromwell.database.slick

import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.tables.{JobStoreEntry, JobStoreResultSimpletonEntry}

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreSlickDatabase extends JobStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  private def addJobStoreEntry(jobStoreEntry: JobStoreEntry, jobStoreSimpletonEntries: Int => Iterable[JobStoreResultSimpletonEntry])
                              (implicit ec: ExecutionContext): DBIO[Unit] = {
    for {
      insertedJobStoreId <- dataAccess.jobStoreAutoInc += jobStoreEntry
      _ <- dataAccess.jobStoreResultSimpletonAutoInc ++= jobStoreSimpletonEntries(insertedJobStoreId)
    } yield ()
  }

  override def addJobStoreEntries(entries: Iterable[(JobStoreEntry, Int => Iterable[JobStoreResultSimpletonEntry])])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = DBIO.sequence(entries map (addJobStoreEntry _).tupled)
    runTransaction(action) map { _ => () }
  }

  override def fetchJobResult(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                    (implicit ec: ExecutionContext): Future[Option[(JobStoreEntry, Iterable[JobStoreResultSimpletonEntry])]] = {

    val action = for {
      jobResults <- dataAccess.jobStoreEntriesByJobStoreKey(workflowUuid, callFqn, index, attempt).result.headOption
      simpletons <- jobResults match {
        case Some(r) => dataAccess.jobStoreResultSimpletonsForJobStoreId(r.jobStoreId.get).result
        case _ => DBIO.successful(Seq.empty)
      }
    } yield jobResults map { j => j -> simpletons }
    
    runTransaction(action)
  }

  override def removeJobResultsForWorkflows(workflowUuids: Iterable[String])(implicit ec: ExecutionContext): Future[Int] = {
    val actions = workflowUuids map { workflowUuid => dataAccess.jobStoreEntryByWorkflowUuid(workflowUuid).delete }
    runTransaction(DBIO.sequence(actions)) map { _.sum }
  }
}
